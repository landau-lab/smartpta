params.model_gpu = "/gpfs/commons/home/jzinno/ug-deepvariant/Ultima_parabricks_4.0.4-1.ultimadeepvar2_V100_noTF32.eng"
params.model_cpu = "/gpfs/commons/home/jzinno/concordanz/deepvariant/model/germline/v1.2_rc2/model.ckpt-710000"
params.monkey_dir = "/gpfs/commons/home/jzinno/concordanz/deepvariant/scripts"

process UGDeepVariantGPU {
    if ("${workflow.stubRun}" == "false") {
        memory '56 GB'
        cpus 10
        queue 'gpu'
        accelerator 1
    }

    tag 'ug-deepvariant'

    publishDir "${params.out}/ug-deepvariant", mode: 'symlink'

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("${bam_file.baseName}.g.vcf.gz"), emit: gvcfs
    path("${bam_file.baseName}.g.vcf.gz.tbi"), emit: gvcf_indices


    script:
    """
    module load gcloud
    module load singularity
    module load htslib/1.18

    export SINGULARITY_DOCKER_USERNAME='_token'
    export SINGULARITY_DOCKER_REGISTRY="gcr.io"
    export SINGULARITY_DOCKER_PASSWORD="\$(gcloud auth print-access-token)"

    singularity run --bind \$(dirname ${params.model_gpu}),\$(dirname ${params.ref}),\$(dirname \$(readlink -f ${bam_file})),\$PWD --nv docker://us.gcr.io/nygc-comp-p-f9e9/clara-parabricks:4.1.0-1.ultimamay \
        pbrun deepvariant \
        --ref ${params.ref} \
        --in-bam \$(readlink -f ${bam_file})  \
        --out-variants \$PWD/${bam_file.baseName}.g.vcf \
        --num-gpus 1 \
        --pb-model-file ${params.model_gpu} \
        --channel-hmer-deletion-quality \
        --channel-hmer-insertion-quality \
        --channel-non-hmer-insertion-quality \
        --aux-fields-to-keep tp,t0 \
        --skip-bq-channel \
        --min-base-quality 5 \
        --dbg-min-base-quality 0 \
        --vsc-min-fraction-indels 0.06 \
        --vsc-min-fraction-snps 0.12 \
        --ws-min-windows-distance 20 \
        --max-read-size-512 \
        --no-channel-insert-size \
        --disable-use-window-selector-model \
        --p-error 0.005 \
        --channel-ins-size \
        --max-ins-size 10 \
        --vsc-min-fraction-hmer-indels 0.12 \
        --consider-strand-bias \
        --vsc-turn-on-non-hmer-ins-proxy-support \
        --gpu-num-per-partition 1 \
        --run-partition \
        --num-cpu-threads-per-stream 5 \
        --num-streams-per-gpu 2 \
        --gvcf

    bgzip -@${task.cpus} ${bam_file.baseName}.g.vcf
    tabix -p vcf ${bam_file.baseName}.g.vcf.gz
    """
  stub:
    """
    touch ${bam_file.baseName}.g.vcf.gz
    touch ${bam_file.baseName}.g.vcf.gz.tbi
    """
}

process UGDeepVariantCPU {
    if ("${workflow.stubRun}" == "false") {
        memory '32 GB'
        cpus 10
    }

    tag 'ug-deepvariant'

    publishDir "${params.out}/ug-deepvariant", mode: 'symlink'

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("${bam_file.baseName}.g.vcf.gz"), emit: gvcfs
    path("${bam_file.baseName}.g.vcf.gz.tbi"), emit: gvcf_indices


    script:
    """
    module load gcloud
    module load singularity
    module load htslib/1.18

    export SINGULARITY_DOCKER_USERNAME='_token'
    export SINGULARITY_DOCKER_REGISTRY="gcr.io"
    export SINGULARITY_DOCKER_PASSWORD="\$(gcloud auth print-access-token)"
    export TMPDIR=\$PWD/tmp

    if [ ! -d "\$TMPDIR" ]; then
        mkdir -p "\$TMPDIR"
    fi

    singularity run --bind \$(dirname ${params.model_cpu}),\$(dirname ${params.ref}),\$(dirname \$(readlink -f ${bam_file})),\$PWD docker://us-central1-docker.pkg.dev/ganymede-331016/ultimagen/deepvariant:ug-1.4.13 \
        ${params.monkey_dir}/run_deepvariant \
        --make_examples_extra_args min_base_quality=5,dbg_min_base_quality=0,vsc_min_fraction_indels=0.06,vsc_min_fraction_hmer_indels=0.12,vsc_min_fraction_snps=0.12,vsc_min_count_snps=2,ws_min_windows_distance=20,min_mapping_quality=5,candidate_min_mapping_quality=5,max_reads_per_partition=1500,aux_fields_to_keep=tp+t0,skip_bq_channel=true,channels=hmer_deletion_quality+hmer_insertion_quality+non_hmer_insertion_quality,add_ins_size_channel=true,max_ins_size=10,p_error=0.005,optimal_coverages=50 \
        --model_type WGS \
        --customized_model ${params.model_cpu} \
        --ref ${params.ref} \
        --reads \$(readlink -f ${bam_file}) \
        --num_shards ${task.cpus} \
        --output_vcf \$PWD/${bam_file.baseName}.vcf.gz \
        --output_gvcf \$PWD/${bam_file.baseName}.g.vcf.gz \
        --intermediate_results_dir \$PWD/tmp \

    rm -rf \$PWD/tmp

    """
  stub:
    """
    touch ${bam_file.baseName}.g.vcf.gz
    touch ${bam_file.baseName}.g.vcf.gz.tbi
    """
}

process ILDeepVariant {
    if ("${workflow.stubRun}" == "false") {
        memory '56 GB'
        cpus 10
        queue 'gpu'
        accelerator 1
    }

    tag 'il-deepvariant'

    publishDir "${params.out}/il-deepvariant", mode: 'symlink'

    input:
    path(bam_file)
    path(bam_index)

    output:
    path("${bam_file.baseName}.g.vcf.gz"), emit: gvcfs
    path("${bam_file.baseName}.g.vcf.gz.tbi"), emit: gvcf_indices


    script:
    """
    module load gcloud
    module load singularity
    module load htslib/1.18

    export SINGULARITY_DOCKER_USERNAME='_token'
    export SINGULARITY_DOCKER_REGISTRY="gcr.io"
    export SINGULARITY_DOCKER_PASSWORD="\$(gcloud auth print-access-token)"

    singularity run --bind \$(dirname ${params.ref}),\$(dirname \$(readlink -f ${bam_file})),\$PWD --nv docker://gcr.io/nygc-comp-p-f9e9/clara-parabricks:4.1.1-1 \
        pbrun deepvariant \
        --ref ${params.ref} \
        --in-bam \$(readlink -f ${bam_file})  \
        --out-variants \$PWD/${bam_file.baseName}.g.vcf \
        --num-gpus 1 \
        --gpu-num-per-partition 1 \
        --run-partition \
        --num-cpu-threads-per-stream 5 \
        --num-streams-per-gpu 2 \
        --gvcf

    bgzip -@${task.cpus} ${bam_file.baseName}.g.vcf
    tabix -p vcf ${bam_file.baseName}.g.vcf.gz
    """
  stub:
    """
    touch ${bam_file.baseName}.g.vcf.gz
    touch ${bam_file.baseName}.g.vcf.gz.tbi
    """
}
