#!/usr/bin/env nextflow

include { FastP } from '../modules/fastp'
include { Trust4 } from '../modules/trust4'

workflow {
    Channel
        .fromPath( params.rna_fastq_table )
        .splitCsv(sep:'\s', header:false)
        .map { row ->
            tuple( row[0], row[1] )
        }
        .set { fastq_ch }
    FastP( fastq_ch  )
    Trust4( FastP.out.trimmed )
}