#!/usr/bin/env nextflow


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

    UGDeepVariant.out.gvcfs
        .map { gvcf -> gvcf.toString() }
        .collectFile(name: params.gvcf_list, newLine: true)
        .set { gvcf_list_ch }

    GLNexus(gvcf_list_ch, params.sample_id)
}