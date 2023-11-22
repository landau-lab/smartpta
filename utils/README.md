# Utils

### A collection of ultilities that hopefully will have utility.

## grip

Pull a symlink target into its source

**NEVER USE THIS WHEN A PIPELINE IS RUNNING**

### Usage

```bash
$ grip <file>
#or
$ find <dir> | grip
#or
$ grip <file> <file> <file> ...
#or
$ grip *
#etc
```

For example, if you'd like to delete a piplines work directory, but you don't want to delete the symlinked output, you can use grip to pull the output from work directory, then delete the work directory.

```bash
$ tree output
output
├── dedup
│   ├── cell1.dedup.bam -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/52/4916914fc26b522f2d3e49cada7e07/cell1.dedup.bam
│   ├── cell1.dedup.bam.bai -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/52/4916914fc26b522f2d3e49cada7e07/cell1.dedup.bam.bai
│   ├── cell2.dedup.bam -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/a3/f44a09ee42d78ddd1208fcbe2344b9/cell2.dedup.bam
│   ├── cell2.dedup.bam.bai -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/a3/f44a09ee42d78ddd1208fcbe2344b9/cell2.dedup.bam.bai
│   ├── cell3.dedup.bam -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/fa/4a0bd210aeb8b62a43fe837acfdf26/cell3.dedup.bam
│   └── cell3.dedup.bam.bai -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/fa/4a0bd210aeb8b62a43fe837acfdf26/cell3.dedup.bam.bai
├── glnexus
│   └── test.glnexus.bcf -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/26/6c4d3e65c7ddbca437c5b8362ca510/test.glnexus.bcf
└── ug-deepvariant
    ├── cell1.dedup.g.vcf.gz -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/2f/82c80fe7b1780a6d4eeb4eccca0ae0/cell1.dedup.g.vcf.gz
    ├── cell1.dedup.g.vcf.gz.tbi -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/2f/82c80fe7b1780a6d4eeb4eccca0ae0/cell1.dedup.g.vcf.gz.tbi
    ├── cell2.dedup.g.vcf.gz -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/20/25b2c14ecc87a70dee26b6eac81b24/cell2.dedup.g.vcf.gz
    ├── cell2.dedup.g.vcf.gz.tbi -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/20/25b2c14ecc87a70dee26b6eac81b24/cell2.dedup.g.vcf.gz.tbi
    ├── cell3.dedup.g.vcf.gz -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/a2/4a9972167202614d1fb8c891ad0c1f/cell3.dedup.g.vcf.gz
    └── cell3.dedup.g.vcf.gz.tbi -> /gpfs/commons/groups/landau_lab/jzinno/test-nf/darkshore/work/a2/4a9972167202614d1fb8c891ad0c1f/cell3.dedup.g.vcf.gz.tbi


$ find output | grip

$ tree output

output
├── dedup
│   ├── cell1.dedup.bam
│   ├── cell1.dedup.bam.bai
│   ├── cell2.dedup.bam
│   ├── cell2.dedup.bam.bai
│   ├── cell3.dedup.bam
│   └── cell3.dedup.bam.bai
├── glnexus
│   └── test.glnexus.bcf
└── ug-deepvariant
    ├── cell1.dedup.g.vcf.gz
    ├── cell1.dedup.g.vcf.gz.tbi
    ├── cell2.dedup.g.vcf.gz
    ├── cell2.dedup.g.vcf.gz.tbi
    ├── cell3.dedup.g.vcf.gz
    └── cell3.dedup.g.vcf.gz.tbi
```
