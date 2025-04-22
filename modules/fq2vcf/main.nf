process FastqToVCF {
    if ("${workflow.stubRun}" == "false") {
        memory '110 GB'
        cpus 8
        accelerator 1
        clusterOptions '--gres gpu:1'
    }

    tag 'fq2vcf'

    container 'docker://zinno/parabricks:4.2.1-1b'

    publishDir "${params.out}/fq2vcf", mode: 'symlink'

    input:
    path(fastqs)

    output:
    path("*.bam"), emit: bam
    path("*.bam.bai"), emit: bam_indices
    path("*.g.vcf.gz"), emit: gvcfs
    path("*.g.vcf.gz.tbi"), emit: gvcf_indices


    script:
    """
    prefix=\$(echo "${fastqs[0].simpleName}" | rev | cut -d'_' -f2- | rev)

    pbrun deepvariant_germline \
        --ref ${params.ref} \
        --in-fq ${fastqs[0]} ${fastqs[1]} \
        --out-bam \${prefix}.bam \
        --out-variants \${prefix}.g.vcf \
        --read-group-sm \${prefix} \
        --num-gpus 1 \
        --gvcf

    bgzip -@${task.cpus} \${prefix}.g.vcf
    tabix -p vcf \${prefix}.g.vcf.gz
    """
  stub:
    """
    prefix=\$(echo "${fastqs[0].simpleName}" | rev | cut -d'_' -f2- | rev)
    touch \${prefix}.bam
    touch \${prefix}.bam.bai
    touch \${prefix}.g.vcf.gz
    touch \${prefix}.g.vcf.gz.tbi
    """
}