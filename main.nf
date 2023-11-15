#!/usr/bin/env nextflow

params.bam_list = "bams.txt"  // A text file with the path to each BAM file on a separate line
params.ref = "path/to/reference.fasta"
params.model = "/gpfs/commons/home/jzinno/ug-deepvariant/Ultima_parabricks_4.0.4-1.ultimadeepvar2_V100_noTF32.eng"
params.gvcf_list = "gvcfs.txt"  // A text file with the path to each GVCF file on a separate line
params.out = "output"
params.sample_id = "test"

include { MarkDuplicatesSpark } from './modules/dedup'
include { UGDeepVariant } from './modules/deepvariant'
include { GLNexus } from './modules/glnexus'

workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }

    MarkDuplicatesSpark(bams_ch.map { it }, params.ref)
    UGDeepVariant(MarkDuplicatesSpark.out.dedup_bam, params.ref)
    GLNexus(UGDeepVariant.out.gvcfs.collectFile(name: params.gvcf_list, newLine: true), params.sample_id)
}