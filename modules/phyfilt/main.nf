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


    script:
    """
    module load bcftools/1.18
    bcftools filter -e 'INFO/avsnp150!="." || F_MISSING > 0.5 || QUAL < 30 || TYPE != "snp" || N_ALT > 1 || AC < 2 || AC > 1000 || AF < 0.01 || AF > 0.99 || STRLEN(REF) != 1 || STRLEN(ALT) != 1' ${annotated_vcf} -Oz -o ${annotated_vcf.simpleName}.phyfilt.vcf.gz
    tabix -p vcf ${annotated_vcf.simpleName}.phyfilt.vcf.gz

    """
    stub:
    """
    touch ${annotated_vcf.simpleName}.phyfilt.vcf.gz
    touch ${annotated_vcf.simpleName}.phyfilt.vcf.gz.tbi
    """

}