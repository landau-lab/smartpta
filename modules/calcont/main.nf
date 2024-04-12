process CalcCont {
    if ("${workflow.stubRun}" == "false") {
        memory '2 GB'
        cpus 1
    }

    tag 'crosscontamination'

    container 'docker://broadinstitute/gatk:4.5.0.0'

    publishDir "${params.out}/crosscont", mode: 'symlink'

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("${bam_file.baseName}.crosscont.table")

    script:
    """
    gatk GetPileupSummaries \
        -I ${bam_file} \
        -V ${params.resource_dir}/somatic-hg38/small_exac_common_3.hg38.vcf.gz \
        -L ${params.resource_dir}/somatic-hg38/small_exac_common_3.hg38.vcf.gz \
        -O Getpileupsummaries.${bam_file.baseName}.table

    gatk CalculateContamination \
        -I Getpileupsummaries.${bam_file.baseName}.table \
        -O ${bam_file.baseName}.crosscont.table
    
    """
    stub:
    """
    touch ${bam_file.baseName}.crosscont.table
    """
}
