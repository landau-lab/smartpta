#!/bin/bash
#SBATCH --job-name cpu-ug-dv
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=56G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jzinno@nygenome.org

module load gcloud
module load singularity
module load htslib/1.18

export SINGULARITY_DOCKER_USERNAME='_token'
export SINGULARITY_DOCKER_REGISTRY="gcr.io"
export SINGULARITY_DOCKER_PASSWORD="$(gcloud auth print-access-token)"

ref=/gpfs/commons/groups/landau_lab/ug-sc1000-dev/ref/hg38/v0/Homo_sapiens_assembly38.fasta
model=/gpfs/commons/home/jzinno/concordanz/deepvariant/model/germline/v1.2_rc2/model.ckpt-710000
out=$PWD

export TMPDIR=$out/tmp

if [ ! -d "$TMPDIR" ]; then
    mkdir -p "$TMPDIR"
fi

sample=$1
identifier=$(basename $sample .bam)

singularity run --bind $(dirname $model),$(dirname $ref),$(dirname $sample),$out docker://us-central1-docker.pkg.dev/ganymede-331016/ultimagen/deepvariant:ug-1.4.13 \
    /opt/deepvariant/bin/run_deepvariant \
    --make_examples_extra_args min_base_quality=5,dbg_min_base_quality=0,vsc_min_fraction_indels=0.06,vsc_min_fraction_hmer_indels=0.12,vsc_min_fraction_snps=0.12,vsc_min_count_snps=2,ws_min_windows_distance=20,min_mapping_quality=5,candidate_min_mapping_quality=5,max_reads_per_partition=1500,aux_fields_to_keep=tp+t0,skip_bq_channel=true,channels=hmer_deletion_quality+hmer_insertion_quality+non_hmer_insertion_quality,add_ins_size_channel=true,max_ins_size=10,p_error=0.005,optimal_coverages=50 \
    --model_type WGS \
    --customized_model $model \
    --ref $ref \
    --reads $sample \
    --num_shards 8 \
    --sample_name $identifier \
    --output_vcf $out \
    --output_gvcf $out \
    --intermediate_results_dir $out/tmp \
