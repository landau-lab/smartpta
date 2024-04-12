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
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz"), emit: annovar_vcf
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi")

    script:
    """
    module load bcftools/1.18

    ${params.annovar_path}/table_annovar.pl \
        ${variant_file} \
        ${params.resource_dir}/humandb/ \
        -buildver hg38 \
        -out ${variant_file.simpleName} \
        -protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
        -remove \
        -operation g,f,f,f,f,f \
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

process AnnovarRAMDisk {
    if ("${workflow.stubRun}" == "false") {
        memory '256 GB'
        cpus 6
        queue 'bigmem'
    }
    tag "annotation"

    publishDir "${params.out}/annovar", mode: 'symlink'

    input:
    path(variant_file)

    output:
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz"), emit: annovar_vcf
    path("${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi")

    script:
    """
    # Setting up the RAM disk
    mount -t tmpfs -o size=72G tmpfs /tmp/humandbRAMdisk/
    cp -r ${params.resource_dir}/humandb/ /tmp/humandbRAMdisk/

    # Running ANNOVAR with the database in the RAM disk
    module load bcftools/1.18

    ${params.annovar_path}/table_annovar.pl \
        ${variant_file} \
        /tmp/humandbRAMdisk/ \
        -buildver hg38 \
        -out ${variant_file.simpleName} \
        -protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
        -remove \
        -operation g,f,f,f,f,f \
        -nastring . \
        -vcfinput \
        -thread ${task.cpus}

    bgzip -@${task.cpus} ${variant_file.simpleName}.hg38_multianno.vcf
    tabix -p vcf ${variant_file.simpleName}.hg38_multianno.vcf.gz
    rm ${variant_file.simpleName}.avinput
    rm ${variant_file.simpleName}.hg38_multianno.txt

    # Cleanup
    rm -r /tmp/humandbRAMdisk/
    umount /tmp/humandbRAMdisk/
    """
    stub:
    """
    touch ${variant_file.simpleName}.hg38_multianno.vcf.gz
    touch ${variant_file.simpleName}.hg38_multianno.vcf.gz.tbi
    """

}

process AnnovarOCI {
    if ("${workflow.stubRun}" == "false") {
        memory '256 GB'
        cpus 6
        queue 'bigmem'
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
        -protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
        -remove \
        -operation g,f,f,f,f,f \
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