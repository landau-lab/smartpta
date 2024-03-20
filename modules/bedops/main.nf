process GenerateIntervals {
    if ("${workflow.stubRun}" == "false") {
        memory '2 GB'
        cpus 1
        queue 'pe2'
    }

    tag 'bedops'

    publishDir "${params.out}/bedops", mode: 'symlink'

    input:
    path(params.ref_idx)

    output:
    path("intervals.bed"), emit: intervals

    script:
    """
    module load bedops 2.4.35

    awk '{print \$1"\t"0"\t"\$2}' ${params.ref_idx} | head -n24 | sort-bed - | bedops --chop ${params.chop} - | sort --version-sort |  awk '{print \$1":"\$2"-"\$3}' > intervals.bed


    """
    stub:
    """
    module load bedops 2.4.35

    awk '{print \$1"\t"0"\t"\$2}' ${params.ref_idx} | head -n24 | sort-bed - | bedops --chop ${params.chop} - | sort --version-sort |  awk '{print \$1":"\$2"-"\$3}' > intervals.bed
    """
}

process SplitVCF {
    if ("${workflow.stubRun}" == "false") {
        memory '2 GB'
        cpus 1
        queue 'pe2'
    }

    tag 'bedops'

    publishDir "${params.out}/chunks", mode: 'symlink'

    input:
    path(joint_vcf)
    val(interval)

    output:
    path("${joint_vcf.baseName}_${interval}.vcf.gz"), emit: split_vcf

    script:
    """
    module load bcftools/1.18
    bcftools view -r ${interval} ${joint_vcf} -Oz > ${joint_vcf.baseName}_${interval}.vcf.gz
    tabix -p vcf ${joint_vcf.baseName}_${interval}.vcf.gz
    """
    stub:
    """
    touch ${joint_vcf.baseName}_${interval}.vcf.gz
    touch ${joint_vcf.baseName}_${interval}.vcf.gz.tbi
    """

}