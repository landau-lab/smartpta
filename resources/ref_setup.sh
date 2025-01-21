#!/bin/bash
set -euo pipefail

log() {
  echo "$(date): $1" 
}

for cmd in gsutil curl; do
  if ! command -v $cmd &> /dev/null; then
    log "Error: $cmd command not found."
    exit 1
  fi
done

log "Downloading hg38 reference bundle..."
gsutil -m cp -r "gs://genomics-public-data/references/hg38" . || { log "Error downloading hg38 bundle"; exit 1; }

log "Downloading somatic best practices resources..."  
gsutil -m cp -r "gs://gatk-best-practices/somatic-hg38" . || { log "Error downloading somatic resources"; exit 1; }

log "Downloading and setting up ANNOVAR databases..."
mkdir -p humandb
cd humandb
declare -a dbs=("avsnp150" "clinvar_20220320" "cosmic70" "dbnsfp42c" "exac03" "gnomad40_genome" "gnomad40_exome")
for db in "${dbs[@]}"; do
  log "Downloading $db database..."
  curl -s "http://www.openbioinformatics.org/annovar/download/hg38_${db}.txt.gz" | gunzip -c > "hg38_${db}.txt" || { log "Error downloading or decompressing $db"; exit 1; }
  log "Downloading index for $db database..."
  curl -s "http://www.openbioinformatics.org/annovar/download/hg38_${db}.txt.idx.gz" | gunzip -c > "hg38_${db}.txt.idx" || { log "Error downloading or decompressing index of $db"; exit 1; }
done
log "Downloading refGene database..."
declare -a fa_dbs=("hg38_refGene.txt" "hg38_refGeneMrna.fa" "hg38_refGeneVersion.txt" "hg38_refGeneWithVerMrna.fa" "hg38_refGeneWithVer.txt")
for db in "${fa_dbs[@]}"; do
  log "Downloading $db..."
  curl -s "http://www.openbioinformatics.org/annovar/download/${db}.gz" | gunzip -c > "$db" || { log "Error downloading or decompressing $db"; exit 1; }
done
cd ..

log "Downloading RNA references from 10X (2020 version)..."
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz" || { log "Error downloading 10X RNA refs (2020)"; exit 1; }
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz || { log "Error extracting 10X RNA refs (2020)"; exit 1; }
rm refdata-gex-GRCh38-2020-A.tar.gz

log "Downloading RNA references from 10X (2024 version)..."
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" || { log "Error downloading 10X RNA refs (2024 version)"; exit 1; }
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz || { log "Error extracting 10X RNA refs (2024 version)"; exit 1; }
gunzip refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz || { log "Error unzipping genes.gtf (2024 version)"; exit 1; }
rm refdata-gex-GRCh38-2024-A.tar.gz


log "Reference setup completed successfully"
