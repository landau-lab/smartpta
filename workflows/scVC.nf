#!/usr/bin/env nextflow

include { FlowMarkDuplicates } from '../modules/dedup'
include { UGDeepVariant } from '../modules/deepvariant'
include { GLNexus } from '../modules/glnexus'
include { SingleCheck } from '../modules/singlecheck'

workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }

    FlowMarkDuplicates(bams_ch.map { it })
    UGDeepVariant(FlowMarkDuplicates.out.dedup_bam, params.ref)
    SingleCheck(FlowMarkDuplicates.out.dedup_bam)
    UGDeepVariant.out.gvcfs
        .map { gvcf -> gvcf.toString() }
        .collectFile(name: params.gvcf_list, newLine: true)
        .set { gvcf_list_ch }

    GLNexus(gvcf_list_ch, params.sample_id)
}