#!/usr/bin/env nextflow

include { SplitAnno } from './subworkflows/splitAnno.nf'

workflow Anno {
    SplitAnno(params.joint_vcf)
}

