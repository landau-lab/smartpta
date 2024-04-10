#!/usr/bin/env nextflow

include { SplitAnno } from './subworkflows/splitAnno.nf'

workflow {
    SplitAnno(Channel.fromPath(params.joint_vcf))
}

