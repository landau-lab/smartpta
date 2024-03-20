#!/usr/bin/env nextflow

include { Annovar } from '../modules/annovar'
include { GenerateIntervals; SplitVCF } from '../modules/bedops'

workflow {
    GenerateIntervals(params.ref_idx)
    GenerateIntervals.out.intervals
        .splitText()
        .set { intervals_chopped }
    SplitVCF( params.joint_vcf, intervals_chopped)
    Annovar(SplitVCF.out.split_vcf)

}
    




