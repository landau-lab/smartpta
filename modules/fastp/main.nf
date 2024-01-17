params.adapters = "${moduleDir}/adapters.fasta"

process FastP {
    if ("${workflow.stubRun}" == "false") {
        memory "8 GB"
        cpus 4
        queue "pe2"
    }
    tag "trimming"

    publishDir "${params.out}/fastp", mode: 'symlink'

    input:
    path(fastqs)


    output:
    tuple path("${fastqs[0].simpleName}.fastp.fastq.gz"), path("${fastqs[1].simpleName}.fastp.fastq.gz"), emit: trimmed
    path("*.fastp.json"), emit: fastp_json
    path("*.fastp.html")
    path("*.fastp.failed.fastq.gz"), emit: failed


    script:
    """
    module load fastp/0.23.1

    prefix=\$(echo "${fastqs[0].simpleName}" | rev | cut -d'_' -f2- | rev)

    fastp \
        -i ${fastqs[0]} \
        -I ${fastqs[1]} \
        -o ${fastqs[0].simpleName}.fastp.fastq.gz \
        -O ${fastqs[1].simpleName}.fastp.fastq.gz \
        --overrepresentation_analysis \
        --adapter_fasta ${params.adapters} \
        --average_qual 20 \
        --json \$prefix.fastp.json \
        --html \$prefix.fastp.html \
        --failed_out \$prefix.fastp.failed.fastq.gz
        --thread ${task.cpus} \

    """
    stub:
    """
    prefix=\$(echo "${fastqs[0].simpleName}" | rev | cut -d'_' -f2- | rev)
    touch ${fastqs[0].simpleName}.fastp.fastq.gz
    touch ${fastqs[1].simpleName}.fastp.fastq.gz
    touch \${prefix}.fastp.json
    touch \${prefix}.fastp.html
    """
}