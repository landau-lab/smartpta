process FastP {
    if ("${workflow.stubRun}" == "false") {
        memory "8 GB"
        cpus 4
        queue "pe2"
    }
    tag "fastp"

    publishDir "${params.out}/fastp", mode: 'symlink'

    input:
    path(fastqs)


    output:
    tuple path("${fastqs[0].simpleName}.fastp.fastq.gz"), path("${fastqs[1].simpleName}.fastp.fastq.gz"), emit: trimmed
    path("${fastqs[0].simpleName}.fastp.json")
    path("${fastqs[0].simpleName}.fastp.html")


    script:
    """
    module load fastp/0.23.1

    fastp \
        -i ${fastqs[0]} \
        -I ${fastqs[1]} \
        -o ${fastqs[0].simpleName}.fastp.fastq.gz \
        -O ${fastqs[1].simpleName}.fastp.fastq.gz \
        --overrepresentation_analysis \
        --adapter_fasta ${params.adapters} \
        --average_qual 20 \
        --json ${fastqs[0].simpleName}.fastp.json \
        --html ${fastqs[0].simpleName}.fastp.html \
        --thread ${task.cpus} \


    """
    stub:
    """
    touch ${fastqs[0].simpleName}.fastp.fastq.gz
    touch ${fastqs[1].simpleName}.fastp.fastq.gz
    touch ${fastqs[0].simpleName}.fastp.json
    touch ${fastqs[0].simpleName}.fastp.html
    """
}