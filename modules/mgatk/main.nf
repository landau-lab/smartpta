process MitoCall {
    if ("${workflow.stubRun}" == "false") {
        memory '5 GB'
        cpus 12
        queue 'pe2'
    }

    conda "${moduleDir}/env.yml"

    tag 'mgatk'

    publishDir "${params.out}/mgatk", mode: 'symlink'

    input:
    path(mgatk_bams)

    output:
    path("final/${params.sample_id}.rds")


    script:
    """
    mgatk call -i \$PWD \
            --mito-genome 'hg38' \
            --output \$PWD \
            --name ${params.sample_id} \
            --ncores ${task.cpus} \
            --alignment-quality 20 \
            --base-qual 20 \
            --emit-base-qualities
    """
    stub:
    """
    mkdir final
    touch final/${params.sample_id}.rds
    """
}