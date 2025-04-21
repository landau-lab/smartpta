#!/usr/bin/env nextflow

include { FastqToVCF } from '../modules/fq2vcf'
include { CalcCont } from '../modules/calcont'
include { GLNexus } from '../modules/glnexus'
include { SplitAnno } from './subworkflows/splitAnno.nf'


workflow {
    Channel
        .fromPath( params.fq_list )
        .splitCsv(sep:'\s', header:false)
        .map { row ->
            tuple( row[0], row[1] )
        }
        .set { fastq_ch }
    FastqToVCF( fastq_ch )
    FastqToVCF.out.gvcfs
        .map { gvcf -> gvcf.toString() }
        .collectFile(name: 'gvcfs.txt', newLine: true)
        .set { gvcf_list_ch }
    GLNexus(gvcf_list_ch)
    SplitAnno(GLNexus.out.joint_vcf)
}