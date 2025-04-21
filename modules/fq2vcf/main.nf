process FastqToVCF {
    if ("${workflow.stubRun}" == "false") {
        memory '56 GB'
        cpus 10
        accelerator 1
        clusterOptions '--gres gpu:1'
    }

    tag 'fq2vcf'

    container 'docker://zinno/parabricks:4.5.0-1b'

    publishDir "${params.out}/fq2vcf", mode: 'symlink'

    input:
    path(fastqs)

    output:
    path("${fastqs[0].baseName}.bam"), emit: bam
    path("${fastqs[0].baseName}.g.vcf.gz"), emit: gvcfs
    path("${fastqs[0].baseName}.g.vcf.gz.tbi"), emit: gvcf_indices


    script:
    """
    pbrun deepvariant_germline \
        --ref ${params.ref} \
        --in-fq ${fastqs[0]} ${fastqs[1]} \
        --out-bam ${fastqs[0].baseName}.bam \
        --out-variants ${fastqs[0].baseName}.g.vcf \
        --num-gpus ${task.accelerator} \
        --gpu-num-per-partition 1 \
        --run-partition \
        --num-cpu-threads-per-stream 5 \
        --num-streams-per-gpu 2 \
        --gvcf

    bgzip -@${task.cpus} ${fastqs[0].baseName}.g.vcf
    tabix -p vcf ${fastqs[0].baseName}.g.vcf.gz
    """
  stub:
    """
    touch ${fastqs[0].baseName}.bam
    touch ${fastqs[0].baseName}.g.vcf.gz
    touch ${fastqs[0].baseName}.g.vcf.gz.tbi
    """
}