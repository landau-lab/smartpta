process GLNexus {
    if ("${workflow.stubRun}" == "false") {
        memory "512 GB"
        cpus 32
    }
    tag "glnexus"

    publishDir "${params.out}/glnexus", mode: 'symlink'

    input:
    path(gvcf_list)


    output:
    path("${params.sample_id}.glnexus.vcf.gz"), emit: joint_vcf
    path("${params.sample_id}.glnexus.vcf.gz.tbi")

    script:
    """
    module load glnexus


    #find unique sample directories
    gvcf_dirs=\$(awk --field-separator '/' 'BEGIN{OFS="/"}{\$NF=""; print \$0}' ${gvcf_list} | sort | uniq | paste -sd,)
    
    cp ${moduleDir}/darkshore.glnexus.yml .

    singularity run --bind \$(dirname \$(readlink -f ${gvcf_list})),\$PWD,\$gvcf_dirs /nfs/sw/glnexus/glnexus-1.4.1/glnexus_cli.sif \
        glnexus_cli \
        --config darkshore.glnexus.yml \
        --list ${gvcf_list} \
        --threads ${task.cpus} \
        --mem-gbytes ${task.memory.toGiga()} \
        > ${params.sample_id}.glnexus.bcf

    module load bcftools/1.18

    bcftools view ${params.sample_id}.glnexus.bcf -Oz > ${params.sample_id}.glnexus.vcf.gz
    tabix -p vcf ${params.sample_id}.glnexus.vcf.gz
    rm ${params.sample_id}.glnexus.bcf
    rm -rf GLnexus.DB

    """
    stub:
    """
    touch ${params.sample_id}.glnexus.vcf.gz
    touch ${params.sample_id}.glnexus.vcf.gz.tbi
    """
}

process GLNexusOCI {
    if ("${workflow.stubRun}" == "false") {
        memory "512 GB"
        cpus 32
    }
    tag "glnexus"

    container 'docker://zinno/glnexus:latest'

    publishDir "${params.out}/glnexus", mode: 'symlink'

    input:
    path(gvcf_list)


    output:
    path("${params.sample_id}.glnexus.vcf.gz"), emit: joint_vcf
    path("${params.sample_id}.glnexus.vcf.gz.tbi")

    script:
    """
    glnexus_cli \
        --config ${moduleDir}/darkshore.glnexus.yml \
        --list ${gvcf_list} \
        --threads ${task.cpus} \
        --mem-gbytes ${task.memory.toGiga()} \
        > ${params.sample_id}.glnexus.bcf


    bcftools view ${params.sample_id}.glnexus.bcf -Oz > ${params.sample_id}.glnexus.vcf.gz
    tabix -p vcf ${params.sample_id}.glnexus.vcf.gz
    rm ${params.sample_id}.glnexus.bcf
    rm -rf GLnexus.DB

    """
    stub:
    """
    touch ${params.sample_id}.glnexus.vcf.gz
    touch ${params.sample_id}.glnexus.vcf.gz.tbi
    """
}