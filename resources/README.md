# Setting up resources

To set up the resources we include a script here that will configure them in the current directory.

```
./ref_setup.sh
```

We use a copy of the `hg38` reference bundle. You can download it via:

```
gsutil -m cp -r "gs://genomics-public-data/references/hg38" .
```

We also use resources from the somatic best practices bucket:

```
gsutil -m cp -r "gs://gatk-best-practices/somatic-hg38" .
```

We also want to download the hg38 annovar dbs:
avsnp150, clinvar_20220320, cosmic70, dbnsfp42c, exac03,and refGene

```
curl -O "http://www.openbioinformatics.org/annovar/download/${db}.txt.gz"
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
├── humandb
├── refdata-gex-GRCh38-2024-A
└── somatic-hg38
```

Once you are done setting up the resources, be sure to update the `nextflow.config` file with the correct paths for the resources. You just need to provide the path to the `ref_bundle` variable at the top of the file.
