params.annovar_path = "/gpfs/commons/groups/landau_lab/tprieto/apps/annovar/"

process Annovar {
    if ("${workflow.stubRun}" == "false") {
        memory '16 GB'
        cpus 4
    }
    tag "annotation"

    container 'docker://zinno/annovar:latest'

    publishDir "${params.out}/annovar", mode: 'symlink'

    input:
    path(variant_file)

    output:
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz"), emit: annovar_vcf
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi")

    script:
    """
    table_annovar.pl \
        ${variant_file} \
        ${params.resource_dir}/humandb/ \
        -buildver hg38 \
        -out ${variant_file.simpleName} \
        -protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320,gnomad40 \
        -remove \
        -operation g,f,f,f,f,f,f \
        -nastring . \
        -vcfinput \
        -thread ${task.cpus}

    bgzip -@${task.cpus} ${variant_file.simpleName}.hg38_multianno.vcf
    tabix -p vcf ${variant_file.simpleName}.hg38_multianno.vcf.gz
    rm ${variant_file.simpleName}.avinput
    rm ${variant_file.simpleName}.hg38_multianno.txt
    """
    stub:
    """
    touch ${variant_file.simpleName}.hg38_multianno.vcf.gz
    touch ${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi
    """

}
