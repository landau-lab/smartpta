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
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.bestModel
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.bestTree
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.bootstraps
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.log
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.mlTrees
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.rba
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.startTree
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.support
    """

}

process CellPhySingleML {
    if ("${workflow.stubRun}" == "false") {
        memory '512 GB'
        cpus 36
        queue 'bigmem'
    }
    tag "phylo"

    publishDir "${params.out}/cellphy/mltrees", mode: 'symlink'

    input:
    path(phylo_vcf)
    each tree_search_idx

    output:
    tuple path("${phylo_vcf.simpleName}.CellPhy.raxml.bestTree.${tree_search_idx}"), path("loglikelihood.${tree_search_idx}.txt")


    script:
    """
    module load cellphy/0.9.2
    raxml-ng-cellphy-linux \
        --search \
        --msa ${phylo_vcf} \
        --model GTGTR4+G+FO \
        --msa-format VCF \
        --threads ${task.cpus} \
        --prefix ${phylo_vcf.simpleName}.CellPhy \
        --tree rand{1} \

    loglikelihood=\$(grep "Final LogLikelihood" ${phylo_vcf.simpleName}.CellPhy.raxml.log | awk '{print $3}')
    echo \$loglikelihood > loglikelihood.${tree_search_idx}.txt

    mv ${phylo_vcf.simpleName}.CellPhy.raxml.bestTree ${phylo_vcf.simpleName}.CellPhy.raxml.bestTree.${tree_search_idx}

    """
    stub:
    """
    touch ${phylo_vcf.simpleName}.CellPhy.raxml.bestTree.${tree_search_idx}
    awk -v seed=\$RANDOM 'BEGIN{srand(seed);print -rand()}' > loglikelihood.${tree_search_idx}.txt
    """

}