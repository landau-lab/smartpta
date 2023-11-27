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
    bcftools filter -e 'INFO/avsnp150!="." || F_MISSING > 0.5 || TYPE != "snp" ||  STRLEN(REF) != 1 || STRLEN(ALT) != 1 || MEDIAN(FMT/DP)<8' ${annotated_vcf} -Oz -o ${annotated_vcf.simpleName}.phyfilt.vcf.gz
    tabix -p vcf ${annotated_vcf.simpleName}.phyfilt.vcf.gz

    """
    stub:
    """
    touch ${annotated_vcf.simpleName}.phyfilt.vcf.gz
    touch ${annotated_vcf.simpleName}.phyfilt.vcf.gz.tbi
    """

}