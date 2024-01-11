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
    path("${trimmed[0].simpleName}.Aligned.sortedByCoord.out.bam") into star_bam
    path("${trimmed[0].simpleName}.Aligned.sortedByCoord.out.bam.bai") into star_bai


    script:
    """
    module load star/2.4.2a
    module load samtools/1.18

    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.star_ref} \
    --sjdbGTFfile ${params.star_gtf} \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${trimmed[0].simpleName}. \
    --readFilesCommand zcat \
    --readFilesIn ${trimmed[0]} ${trimmed[1]} \

    samtools index -@ ${task.cpus} ${trimmed[0].simpleName}.Aligned.sortedByCoord.out.bam

    """
    stub:
    """
    touch ${trimmed[0].simpleName}.Aligned.sortedByCoord.out.bam
    touch ${trimmed[0].simpleName}.Aligned.sortedByCoord.out.bam.bai
    """
}