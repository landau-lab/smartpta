# Setting up resources

We use a copy of the `hg38` reference bundle. You can download it via:

```
gsutil -m cp -r "gs://genomics-public-data/references/hg38" .
```

We also use resources from the somatic best practices bucket:

```
gsutil -m cp -r "gs://gatk-best-practices/somatic-hg38" .
```

To included references for RNA workflows we can pull them from 10X:

```
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
gunzip refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
rm refdata-gex-GRCh38-2024-A.tar.gz
```

The resource dir then should look like:

```
.
├── hg38
├── refdata-gex-GRCh38-2024-A
└── somatic-hg38
```
