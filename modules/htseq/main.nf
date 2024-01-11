process HTSeq {
    if ("${workflow.stubRun}" == "false") {
        memory "4 GB"
        cpus 1
        queue "pe2"
    }
    tag "count"

    publishDir "${params.out}/htseq", mode: 'symlink'

    input:
    path(bam)

    output:
    path("${bam.simpleName}.counts")

    script:
    """
    module load HTSeq/0.6.1

    htseq-count \
        --format bam \
        --order pos \
        --stranded yes \
        --minaqual 10 \
        --type exon \
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