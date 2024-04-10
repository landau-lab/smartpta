process PreFilter {
    if ("${workflow.stubRun}" == "false") {
        memory '16 GB'
        cpus 4
    }
    tag "filtering"

    publishDir "${params.out}/prefilter", mode: 'symlink'

    input:
    path(annotated_vcf)

    output:
    path("${annotated_vcf.simpleName}.prefilter.vcf.gz"), emit: prefilter_vcf
    path("${annotated_vcf.simpleName}.prefilter.vcf.gz.tbi"), emit: prefilter_index


    script:
    """
    module load bcftools/1.18
    bcftools filter -e 'QUAL < 30 || F_MISSING > 0.5 || TYPE != "snp" ||  STRLEN(REF) != 1 || STRLEN(ALT) != 1  || MEDIAN(FMT/DP)<8' ${annotated_vcf} -Oz -o ${annotated_vcf.simpleName}.prefilter.vcf.gz
    tabix -p vcf ${annotated_vcf.simpleName}.prefilter.vcf.gz

    """
    stub:
    """
    touch ${annotated_vcf.simpleName}.prefilter.vcf.gz
    touch ${annotated_vcf.simpleName}.prefilter.vcf.gz.tbi
    """

}

process AnnoFilter {
    if ("${workflow.stubRun}" == "false") {
        memory '16 GB'
        cpus 4
    }
    tag "filtering"

    publishDir "${params.out}/annofilter", mode: 'symlink'

    input:
    path(annotated_vcf)

    output:
    path("${annotated_vcf.simpleName}.annofilter.vcf.gz"), emit: annofilter_vcf
    path("${annotated_vcf.simpleName}.annofilter.vcf.gz.tbi"), emit: annofilter_index


    script:
    """
    module load bcftools/1.18
    bcftools filter -e 'INFO/avsnp150!="."' ${annotated_vcf} -Oz -o ${annotated_vcf.simpleName}.annofilter.vcf.gz
    tabix -p vcf ${annotated_vcf.simpleName}.annofilter.vcf.gz

    """
    stub:
    """
    touch ${annotated_vcf.simpleName}.annofilter.vcf.gz
    touch ${annotated_vcf.simpleName}.annofilter.vcf.gz.tbi
    """

}

process Phyfilt {
    if ("${workflow.stubRun}" == "false") {
        memory '16 GB'
        cpus 4
    }
    tag "filtering"

    publishDir "${params.out}/phyfilt", mode: 'symlink'

    input:
    path(annotated_vcf)

    output:
    path("${annotated_vcf.simpleName}.phyfilt.vcf.gz"), emit: phyfilt_vcf
    path("${annotated_vcf.simpleName}.phyfilt.vcf.gz.tbi"), emit: phyfilt_index


    script:
    """
    module load bcftools/1.18
    bcftools filter -e 'INFO/avsnp150!="." || QUAL < 30 || F_MISSING > 0.5 || TYPE != "snp" ||  STRLEN(REF) != 1 || STRLEN(ALT) != 1  || AC < 2 || MEDIAN(FMT/DP)<8' ${annotated_vcf} -Oz -o ${annotated_vcf.simpleName}.phyfilt.vcf.gz
    tabix -p vcf ${annotated_vcf.simpleName}.phyfilt.vcf.gz

    """
    stub:
    """
    touch ${annotated_vcf.simpleName}.phyfilt.vcf.gz
    touch ${annotated_vcf.simpleName}.phyfilt.vcf.gz.tbi
    """

}