#!/usr/bin/env nextflow

include { FastP } from '../modules/fastp'
include { StarTE } from '../modules/star'
include { RNAMultiQC } from '../modules/multiqc'

workflow {
    Channel
        .fromPath( params.rna_fastq_table )
        .splitCsv(sep:'\s', header:false)
        .map { row ->
            tuple( row[0], row[1] )
        }
        .set { fastq_ch }
    FastP( fastq_ch  )
    StarTE( FastP.out.trimmed )
    FastP.out.fastp_json
        .map { fastp_data -> fastp_data.toString() }
        .collectFile( name: 'fastp_data.txt', newLine: true )
        .set { fastp_data_ch }
    StarTE.out.star_log
        .map { star_log -> star_log.toString() }
        .collectFile( name: 'star_log.txt', newLine: true )
        .set { star_log_ch }
    RNAMultiQC( fastp_data_ch, star_log_ch, )
}