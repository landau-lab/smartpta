# darkshore

[![Docker Build](https://github.com/jzinno/darkshore/actions/workflows/docker-build.yml/badge.svg)](https://github.com/jzinno/darkshore/actions/workflows/docker-build.yml)

### Nextflow pipelines for the worlds biggest single-cell multiomic phylogenies

## Quick start & Test run

```bash
   git clone https://github.com/jzinno/darkshore.git

   cd darkshore

   nextflow workflows/scVC.nf -stub-run -profile stub

   #explore example output
   tree -C output
```

## Installation

### Requirements

These workflows were developed using the following:

- Nextflow 22.10.4+
- Singularity 3.8.6+

The pipelines will automatically pull the required containers when run. Our reference data bundle is required in order to run the pipelines, in order to download the required files run:

```bash
   cd resources
   ./ref_setup.sh
   cd ..
```

More information on the reference data bundle can be found [here](https://github.com/jzinno/darkshore/tree/main/resources).

## Usage

### scVariantCalling (scVC.nf)

This will run the following steps:

- Duplicate Marking
- Contamination Estimate
- UG DeepVariant
- GLNexus joint genotyping
- Variant Annotation
- VCF Filtering

Create a bam list

```bash
   #e.g. bam_list.txt
   /path/to/bam1.bam
   /path/to/bam2.bam
   /path/to/bam3.bam
```

```bash
   nextflow workflows/scVC.nf --bam_list <bam_list> --sample_id <sample_id>
```

### scRNAseq Analysis (scRNA.nf)

This will run the following steps:

- Quality control with FastP
- Alignment with STAR
- Quantification with HTSeq
- Merging of counts
- Quality control report with MultiQC

Create an RNA-seq fastq list

```bash
   #e.g. rna_fastq_pairs.txt
   /path/to/fastq1_R1.fastq.gz /path/to/fastq1_R2.fastq.gz
   /path/to/fastq2_R1.fastq.gz /path/to/fastq2_R2.fastq.gz
   /path/to/fastq3_R1.fastq.gz /path/to/fastq3_R2.fastq.gz
```

```bash
nextflow workflows/scRNA.nf --rna_fastq_table <rna_fastq_pairs> --sample_id <sample_id>
```
