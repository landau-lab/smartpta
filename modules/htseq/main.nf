process HTSeq {
    if ("${workflow.stubRun}" == "false") {
        memory "4 GB"
        cpus 1
    }
    tag "count"

    container 'docker://zinno/rnatools:latest'

    publishDir "${params.out}/htseq", mode: 'symlink'

    input:
    path(bam)

    output:
    path("${bam.simpleName}.counts"), emit: htseq_counts

    script:
    """
    htseq-count \
        --format bam \
        --order pos \
        --stranded yes \
        --minaqual 10 \
        --type ${params.htseq_type} \
        --idattr gene_id \
        --mode union \
        ${bam} \
        ${params.gtf} \
        > ${bam.simpleName}.counts

    """
    stub:
    """
    touch ${bam.simpleName}.counts
    """
}