#!/usr/bin/env nextflow

include { FastP } from '../modules/fastp'
include { Star } from '../modules/star'
include { HTSeq } from '../modules/htseq'
include { MergeCounts } from '../modules/bash'
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
    Star( FastP.out.trimmed )
    HTSeq( Star.out.star_bam )
    HTSeq.out.htseq_counts
        .map { counts -> counts.toString() }
        .collectFile( name: 'htseq_counts.txt', newLine: true )
        .set { htseq_counts_ch }
    MergeCounts( htseq_counts_ch )
    FastP.out.fastp_json
        .map { fastp_data -> fastp_data.toString() }
        .collectFile( name: 'fastp_data.txt', newLine: true )
        .set { fastp_data_ch }
    Star.out.star_log
        .map { star_log -> star_log.toString() }
        .collectFile( name: 'star_log.txt', newLine: true )
        .set { star_log_ch }
    RNAMultiQC( fastp_data_ch, star_log_ch, )
}