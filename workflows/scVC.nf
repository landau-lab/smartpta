#!/usr/bin/env nextflow

include { FlowMarkDuplicates } from '../modules/dedup'
include { UGDeepVariantGPU; UGDeepVariantCPU } from '../modules/deepvariant'
include { CalcCont } from '../modules/calcont'
include { GLNexus } from '../modules/glnexus'
include { PreFilter; AnnoFilter } from '../modules/phyfilt'
include { SplitAnno } from 'subworkflows/splitAnno.nf'


workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }
    FlowMarkDuplicates(bams_ch.map { it })
    CalcCont(FlowMarkDuplicates.out.dedup_bam, FlowMarkDuplicates.out.dedup_bam_index)
    if (params.use_gpu){
        UGDeepVariantGPU(FlowMarkDuplicates.out.dedup_bam, FlowMarkDuplicates.out.dedup_bam_index)
        UGDeepVariantGPU.out.gvcfs
            .map { gvcf -> gvcf.toString() }
            .collectFile(name: 'gvcfs.txt', newLine: true)
            .set { gvcf_list_ch }

    }else{
        UGDeepVariantCPU(FlowMarkDuplicates.out.dedup_bam, FlowMarkDuplicates.out.dedup_bam_index)
        UGDeepVariantCPU.out.gvcfs
            .map { gvcf -> gvcf.toString() }
            .collectFile(name: 'gvcfs.txt', newLine: true)
            .set { gvcf_list_ch }
    }
    GLNexus(gvcf_list_ch)
    PreFilter(GLNexus.out.joint_vcf)
    SplitAnno(PreFilter.out.prefilter_vcf)
    AnnoFilter(SplitAnno.out)
}