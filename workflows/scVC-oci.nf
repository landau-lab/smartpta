#!/usr/bin/env nextflow

include { FlowMarkDuplicates } from '../modules/dedup'
include { UGDeepVariantPB } from '../modules/deepvariant'
include { CalcCont } from '../modules/calcont'
include { GLNexusOCI } from '../modules/glnexus'
include { Annovar } from '../modules/annovar'


workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }
    FlowMarkDuplicates(bams_ch.map { it })
    CalcCont(FlowMarkDuplicates.out.dedup_bam, FlowMarkDuplicates.out.dedup_bam_index)
    UGDeepVariantPB(FlowMarkDuplicates.out.dedup_bam, FlowMarkDuplicates.out.dedup_bam_index)
    UGDeepVariantPB.out.gvcfs
        .map { gvcf -> gvcf.toString() }
        .collectFile(name: 'gvcfs.txt', newLine: true)
        .set { gvcf_list_ch }
    GLNexusOCI(gvcf_list_ch)
    Annovar(GLNexusOCI.out.joint_vcf)
}