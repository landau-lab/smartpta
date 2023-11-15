process MarkDuplicatesSpark {
    memory '64 GB'
    cpus 4
    queue 'pe2'
    tag 'dedup'

    publishDir "${params.out}/dedup", mode: 'symlink'

    input:
    path(bam_file)
    val(reference)

    output:
    path("${bam_file.baseName}.dedup.bam"), emit: dedup_bam
    path("${bam_file.baseName}.dedup.bam.bai"), emit: dedup_bam_index

    script:
    """
    module purge
    module load gatk/4.3.0.0
    module unload java
    module load java/1.8
    module load samtools
    
    if [ ! -d $PWD/tmp ]; then
        mkdir $PWD/tmp
    fi

    gatk MarkDuplicatesSpark \
        --flowbased \
        --reference ${reference} \
        --input ${bam_file} \
        --output ${bam_file.baseName}.dedup.bam \
        --conf 'spark.executor.cores=${task.cpus}' \
        --tmp-dir $PWD/tmp \

    samtools index ${bam_file.baseName}.dedup.bam
    """
}