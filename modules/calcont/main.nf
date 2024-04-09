params.res_dir="/gpfs/commons/groups/landau_lab/tprieto/gatk-bundle/hg38"
params.exac_file="bestpractices/small_exac_common_3.hg38.vcf.gz"

process CalcCont {
    if ("${workflow.stubRun}" == "false") {
        memory '2 GB'
        cpus 1
    }

    tag 'crosscontamination'

    publishDir "${params.out}/crosscont", mode: 'symlink'

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("${bam_file.baseName}.crosscont.table")

    script:
    """
    module load gatk/4.1.8.1
    gatk GetPileupSummaries \
        -I ${bam_file} \
        -V ${params.res_dir}/${params.exac_file} \
        -L ${params.res_dir}/${params.exac_file} \
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

process CalcContOCI {
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
        -V ${params.res_dir}/${params.exac_file} \
        -L ${params.res_dir}/${params.exac_file} \
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