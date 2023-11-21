#!/usr/bin/env nextflow

include { FlowMarkDuplicates } from '../modules/dedup'
include { UGDeepVariant } from '../modules/deepvariant'
include { SingleCheck } from '../modules/singlecheck'
include { GLNexus } from '../modules/glnexus'
include { Annovar } from '../modules/annovar'
include { Phyfilt } from '../modules/phyfilt'
include { MitoCall } from '../modules/mgatk'
include { CellPhy } from '../modules/phylo'

workflow {
    Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { row -> file(row.trim()) }
        .set { bams_ch }

    FlowMarkDuplicates(bams_ch.map { it })
    if (params.call_mito){
        MitoCall(FlowMarkDuplicates.out.dedup_bam.collect())
    }
    UGDeepVariant(FlowMarkDuplicates.out.dedup_bam, params.ref)
    SingleCheck(FlowMarkDuplicates.out.dedup_bam, FlowMarkDuplicates.out.dedup_bam_index)
    UGDeepVariant.out.gvcfs
        .map { gvcf -> gvcf.toString() }
        .collectFile(name: params.gvcf_list, newLine: true)
        .set { gvcf_list_ch }

    GLNexus(gvcf_list_ch, params.sample_id)
    Annovar(GLNexus.out.joint_vcf)
    Phyfilt(Annovar.out.annovar_vcf)
    if (params.run_phylo){
        CellPhy(Phyfilt.out.phyfilt_vcf)
    }
}