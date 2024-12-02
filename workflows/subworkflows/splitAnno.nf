#!/usr/bin/env nextflow

include { Annovar } from '../../modules/annovar'
include { GenerateIntervals; SplitVCF; MergeVCFs } from '../../modules/bedops'

workflow SplitAnno {
    take:
    vcf_ch

    main:
    GenerateIntervals(params.ref_idx)
    GenerateIntervals.out.intervals
        .splitText()
        .set { intervals_chopped }
    intervals_chopped
        .combine(vcf_ch)
        .set { intervals_with_vcf }
    SplitVCF(intervals_with_vcf)
    Annovar(SplitVCF.out.split_vcf)
    Annovar.out.annovar_vcf
        .map { vcf -> vcf.toString()}
        .collectFile(name: params.sample_id,  newLine: true)
        .set { annos }
    MergeVCFs(annos)

    emit:
    merged_vcf = MergeVCFs.out.merged_vcf
    merged_vcf_tbi = MergeVCFs.out.merged_vcf_tbi
}
