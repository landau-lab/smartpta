process SingleCheck {
    if ("${workflow.stubRun}" == "false") {
        memory '6 GB'
        cpus 1
        queue 'pe2'
    }

    tag 'dispersion'

    publishDir "${params.out}/singlecheck", mode: 'symlink'

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("${bam_file.baseName}.SingleCheck.txt")

    script:
    """
    export PATH=/gpfs/commons/groups/landau_lab/tprieto/apps/SingleCheck:$PATH

    SingleCheck \
        -N \
        ${bam_file} \
    """
    stub:
    """
    touch ${bam_file.baseName}.SingleCheck.txt
    """
}