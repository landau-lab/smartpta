params.annovar_path = "/gpfs/commons/groups/landau_lab/tprieto/apps/annovar/"

process Annovar {
    if ("${workflow.stubRun}" == "false") {
        memory '16 GB'
        cpus 4
    }
    tag "annotation"

    publishDir "${params.out}/annovar", mode: 'symlink'

    input:
    path(variant_file)

    output:
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz")
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi")

    script:
    """
    module load bcftools/1.18

    ${params.annovar_path}/table_annovar.pl \
        ${variant_file} \
        ${params.annovar_path}/humandb/ \
        -buildver hg38 \
        -out ${variant_file.simpleName} \
        -protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
        -remove \
        -operation g,f,f,f,f,f \
        -nastring . \
        -vcfinput \

    bgzip -@${task.cpus} ${variant_file.simpleName}.hg38_multianno.vcf
    tabix -p vcf ${variant_file.simpleName}.hg38_multianno.vcf.gz
    rm ${variant_file.simpleName}.hg38_multianno.avinput
    rm ${variant_file.simpleName}.hg38_multianno.txt
    """
    stub:
    """
    touch ${variant_file.simpleName}.hg38_multianno.vcf.gz
    touch ${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi
    """

}