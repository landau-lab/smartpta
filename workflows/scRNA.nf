#!/usr/bin/env nextflow

include { FastP } from '../modules/fastp'
include { Star } from '../modules/star'

workflow {
    Channel
        .fromPath( params.rna_fastq_table )
        .splitCsv(sep:'\s', header:false)
        .map { row ->
            tuple( row[0], row[1] )
        }
        .set { fastq_ch }
    FastP( fastq_ch  )
    Star( FastP.out.trimmed )

}