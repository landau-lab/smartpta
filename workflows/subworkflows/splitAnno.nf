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
    interval_chopped
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
    MergeVCFs.out
}
