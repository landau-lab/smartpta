params.star_ref = "/gpfs/commons/groups/landau_lab/rraviram/hg38/refdata-gex-GRCh38-2020-A/star/"
params.star_gtf = "/gpfs/commons/groups/landau_lab/rraviram/hg38/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

process Star {
    if ("${workflow.stubRun}" == "false") {
        memory "32 GB"
        cpus 6
        queue "pe2"
    }
    tag "align"

    publishDir "${params.out}/star", mode: 'symlink'

    input:
    path(trimmed)


    output:
    path("*.Aligned.sortedByCoord.out.bam"), emit: star_bam
    path("*.Aligned.sortedByCoord.out.bam.bai"), emit: star_bai


    script:
    """
    module load star/2.4.2a
    module load samtools/1.18

    prefix=\$(echo "${trimmed[0].simpleName}" | rev | cut -d'_' -f2- | rev)

    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.star_ref} \
    --sjdbGTFfile ${params.star_gtf} \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix \${prefix}. \
    --readFilesCommand zcat \
    --readFilesIn ${trimmed[0]} ${trimmed[1]} \

    samtools index -@ ${task.cpus} \${prefix}.Aligned.sortedByCoord.out.bam

    """
    stub:
    """
    prefix=\$(echo "${trimmed[0].simpleName}" | rev | cut -d'_' -f2- | rev)
    touch \${prefix}.Aligned.sortedByCoord.out.bam
    touch \${prefix}.Aligned.sortedByCoord.out.bam.bai
    """
}