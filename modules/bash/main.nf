process MergeCounts {
    if ("${workflow.stubRun}" == "false") {
        memory "4 GB"
        cpus 1
        queue "pe2"
    }
    tag "collate"

    publishDir "${params.out}/matrix", mode: 'symlink'

    input:
    path(count_list)

    output:
    path("${params.sample_id}.counts.tab")

    script:
    """
    ${moduleDir}/merge_htseq.sh ${count_list} ${params.sample_id}.counts.tab
    """
    stub:
    """
    touch ${params.sample_id}.counts.tab
    """
}