#!/usr/bin/env nextflow

include { FastP } from '../modules/fastp'
include { Star } from '../modules/star'
include { HTSeq } from '../modules/htseq'
include { MergeCounts } from '../modules/bash'

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

}