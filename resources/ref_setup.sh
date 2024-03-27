#!/bin/bash
set -euo pipefail

log() {
  echo "$(date): $1" 
}

# Check if gsutil is available
if ! command -v gsutil &> /dev/null
then
    log "Error: gsutil command not found. Please install Google Cloud SDK and ensure gsutil is in the PATH."
    exit 1
fi

log "Downloading hg38 reference bundle..."
gsutil -m cp -r "gs://genomics-public-data/references/hg38" . || { log "Error downloading hg38 bundle"; exit 1; }

log "Downloading somatic best practices resources..."  
gsutil -m cp -r "gs://gatk-best-practices/somatic-hg38" . || { log "Error downloading somatic resources"; exit 1; }

log "Downloading RNA references from 10X..."
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" || { log "Error downloading 10X RNA refs"; exit 1; }
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz || { log "Error extracting 10X RNA refs"; exit 1; }
gunzip refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz || { log "Error unzipping genes.gtf"; exit 1; }
rm refdata-gex-GRCh38-2024-A.tar.gz

log "Reference setup completed successfully"