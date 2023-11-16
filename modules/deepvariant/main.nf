params.model = "/gpfs/commons/home/jzinno/ug-deepvariant/Ultima_parabricks_4.0.4-1.ultimadeepvar2_V100_noTF32.eng"

process UGDeepVariant {
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
    val(ref)

    output:
    path("${bam_file.baseName}.g.vcf.gz"), emit: gvcfs
    path("${bam_file.baseName}.g.vcf.gz.tbi"), emit: gvcf_indices


    script:
    """
    module load gcloud
    module load singularity
    moodule load htslib/1.18

    export SINGULARITY_DOCKER_USERNAME='_token'
    export SINGULARITY_DOCKER_REGISTRY="gcr.io"
    export SINGULARITY_DOCKER_PASSWORD="\$(gcloud auth print-access-token)"

    singularity run --bind \$(dirname ${params.model}),\$(dirname ${ref}),\$(dirname \$(readlink ${bam_file})) --nv docker://us.gcr.io/nygc-comp-p-f9e9/clara-parabricks:4.1.0-1.ultimamay \
        pbrun deepvariant \
        --ref ${ref} \
        --in-bam \$(readlink ${bam_file})  \
        --out-variants ${bam_file.baseName}.g.vcf \
        --num-gpus 1 \
        --pb-model-file ${params.model} \
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
        --num-cpu-threads-per-stream 2 \
        --num-streams-per-gpu 5 \
        --gvcf

    bgzip ${bam_file.baseName}.g.vcf
    tabix -p vcf ${bam_file.baseName}.g.vcf.gz
    """
  stub:
    """
    touch ${bam_file.baseName}.g.vcf.gz
    touch ${bam_file.baseName}.g.vcf.gz.tbi
    """
}
