#!/usr/bin/env nextflow

include { FlowMarkDuplicatesOCI } from '../modules/dedup'
include { UGDeepVariantPB } from '../modules/deepvariant'
include { CalcContOCI } from '../modules/calcont'
include { GLNexusOCI } from '../modules/glnexus'
include { AnnovarOCI } from '../modules/annovar'


workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }
    FlowMarkDuplicatesOCI(bams_ch.map { it })
    CalcContOCI(FlowMarkDuplicatesOCI.out.dedup_bam, FlowMarkDuplicatesOCI.out.dedup_bam_index)
    UGDeepVariantPB(FlowMarkDuplicatesOCI.out.dedup_bam, FlowMarkDuplicatesOCI.out.dedup_bam_index)
    UGDeepVariantPB.out.gvcfs
        .map { gvcf -> gvcf.toString() }
        .collectFile(name: 'gvcfs.txt', newLine: true)
        .set { gvcf_list_ch }
    GLNexusOCI(gvcf_list_ch)
    AnnovarOCI(GLNexusOCI.out.joint_vcf)
}