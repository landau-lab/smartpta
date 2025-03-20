#!/usr/bin/env nextflow

include { GLNexus } from '../modules/glnexus'
include { SplitAnno } from './subworkflows/splitAnno.nf'

workflow {
    // Create channel from the gvcf list parameter
    Channel
        .fromPath(params.gvcf_list)
        .set { gvcf_list_ch }

    // Run GLNexus
    GLNexus(gvcf_list_ch)

    // Run SplitAnno on GLNexus output
    SplitAnno(GLNexus.out.joint_vcf)
}