#!/usr/bin/env nextflow

include { FlowMarkDuplicates } from '../modules/dedup'
include { MitoCall } from '../modules/mgatk'

workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }
    FlowMarkDuplicates(bams_ch.map { it })
    MitoCall(FlowMarkDuplicates.out.dedup_bam.collect(), FlowMarkDuplicates.out.dedup_bam_index.collect())
}