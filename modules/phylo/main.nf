process CellPhy {
    if ("${workflow.stubRun}" == "false") {
        memory '512 GB'
        cpus 36
        queue 'bigmem'
    }
    tag "phylo"

    publishDir "${params.out}/cellphy", mode: 'symlink'

    input:
    path(phylo_vcf)

    output:
    path("${phylo_vcf.simpleName}.CellPhy.*")


    script:
    """
    module load cellphy/0.9.2
    raxml-ng-cellphy-linux \
        --all \
        --msa ${phylo_vcf} \
        --model GTGTR4+G+FO \
        --msa-format VCF \
        --threads ${task.cpus} \
        --prefix ${phylo_vcf.simpleName}.CellPhy \
        --bs-trees 10 \
        --tree pars{10} \

    """
    stub:
    """
    ${phylo_vcf.simpleName}.CellPhy.raxml.bestModel
    ${phylo_vcf.simpleName}.CellPhy.raxml.bestTree
    ${phylo_vcf.simpleName}.CellPhy.raxml.bootstraps
    ${phylo_vcf.simpleName}.CellPhy.raxml.log
    ${phylo_vcf.simpleName}.CellPhy.raxml.mlTrees
    ${phylo_vcf.simpleName}.CellPhy.raxml.rba
    ${phylo_vcf.simpleName}.CellPhy.raxml.startTree
    ${phylo_vcf.simpleName}.CellPhy.raxml.support
    """

}