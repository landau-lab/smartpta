# darkshore

### Nextflow pipelines for the worlds biggest single-cell multiomic phylogenies

## Setup & Test run

```bash
   git clone https://github.com/jzinno/darkshore.git

   cd darkshore

   module load nextflow/22.10.4

   nextflow workflows/scVC.nf -stub-run -profile stub

   #explore example output
   tree -C
```

Right now only runs on NYGC cluster

## Usage

### scVariantCalling (only pipeline currently available)

This will run the following steps:

- MarkDuplicates
- UG DeepVariant
- GLNexus joint genotyping

Create a bam list

```bash
   #bam_list.txt
   /path/to/bam1.bam
   /path/to/bam2.bam
   /path/to/bam3.bam
```

```bash
   nextflow workflows/scVC.nf --bam_list <bam_list> --sample_id <sample_id>
```
