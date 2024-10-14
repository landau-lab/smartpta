process Trust4 {
    if ("${workflow.stubRun}" == "false") {
        memory "32 GB"
        cpus 6
    }
    tag "align"

    container 'docker://quay.io/biocontainers/trust4:1.0.13--h43eeafb_0'

    publishDir "${params.out}/trust4", mode: 'symlink'

    input:
    path(trimmed)


    output:
    path("*_airr_align.tsv")
    path("*_airr.tsv")
    path("*_annot.fa")
    path("*_assembled_reads.fa")
    path("*_cdr3.out")
    path("*_final.out")
    path("*_raw.out")
    path("*_report.tsv")
    path("*_toassemble_1.fq")
    path("*_toassemble_2.fq")


    script:
    """
    prefix=\$(echo "${trimmed[0].simpleName}" | rev | cut -d'_' -f2- | rev)

    run-trust4 \
    -t ${task.cpus} \
    -1 ${trimmed[0]} \
    -2 ${trimmed[1]} \
    -f ${params.resource_dir}/TRUST4/hg38_bcrtcr.fa \
    --ref ${params.resource_dir}/TRUST4/human_IMGT+C.fa \
    -o \${prefix} \


    """
    stub:
    """
    prefix=\$(echo "${trimmed[0].simpleName}" | rev | cut -d'_' -f2- | rev)
    touch \${prefix}_airr_align.tsv
    touch \${prefix}_airr.tsv
    touch \${prefix}_annot.fa
    touch \${prefix}_assembled_reads.fa
    touch \${prefix}_cdr3.out
    touch \${prefix}_final.out
    touch \${prefix}_raw.out
    touch \${prefix}_report.tsv
    touch \${prefix}_toassemble_1.fq
    touch \${prefix}_toassemble_2.fq
    """
}