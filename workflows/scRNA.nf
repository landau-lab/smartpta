#!/usr/bin/env nextflow

include { FastP } from '../modules/fastp'

workflow {
    Channel
        .fromPath( params.rna_fastq_table )
        .splitCsv()
        .map { row ->
            tuple( row.fastq_1, row.fastq_2 )
        }
        .set { fastq_ch }
    FastP( fastq_ch  )
}