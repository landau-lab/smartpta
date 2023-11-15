process GLNexus {
    if ("${workflow.stubRun}" == "false") {
        memory "512 GB"
        cpus 32
        queue "bigmem"
    }
    tag "glnexus"

    publishDir "${params.out}/glnexus", mode: 'symlink'

    input:
    path(gvcf_list)
    val sample_id

    output:
    path("${sample_id}.glnexus.bcf")

    script:
    """
    module load glnexus


    #find unique sample directories
    gvcf_dirs=\$(awk --field-separator '/' 'BEGIN{OFS="/"}{\$NF=""; print \$0}' ${gvcf_list} | sort | uniq | paste -sd,)

    singularity run --bind \$(dirname ${gvcf_list}),\$gvcf_dirs /nfs/sw/glnexus/glnexus-1.4.1/glnexus_cli.sif \
        glnexus_cli \
        --config DeepVariant \
        --list ${gvcf_list} \
        --bed /gpfs/commons/groups/landau_lab/jzinno/ultima/eso-ug/glnexus/ug_hcr_nohead.bed \
        --threads ${task.cpus} \
        --mem-gbytes ${task.memory.toGiga()} \
        > ${sample_id}.glnexus.bcf
    """
    stub:
    """
    touch ${sample_id}.glnexus.bcf
    """
}