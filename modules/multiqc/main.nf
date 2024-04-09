process RNAMultiQC {
    if ("${workflow.stubRun}" == "false") {
        memory "16 GB"
        cpus 4
    }
    tag "report"

    publishDir "${params.out}/multiqc", mode: 'symlink'

    input:
    path(fastp_data)
    path(star_data)


    output:
    path("${params.sample_id}_multiqc_report.html")

    script:
    """
    module load multiqc/1.19

    cat ${fastp_data} ${star_data} > all_data.txt

    multiqc \
        --file-list all_data.txt \
        --filename ${params.sample_id}_multiqc_report.html


    """
    stub:
    """
    touch ${params.sample_id}_multiqc_report.html
    """
}