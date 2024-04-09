process CellPhy {
    if ("${workflow.stubRun}" == "false") {
        memory '512 GB'
        cpus 36
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

process MLSearchCellPhy {
    if ("${workflow.stubRun}" == "false") {
        memory '4 GB'
        cpus 4
    }
    tag "tree-search"

    publishDir "${params.out}/cellphy/mltrees", mode: 'symlink'

    input:
    path(phylo_vcf)
    each tree_search_idx

    output:
    tuple path("${phylo_vcf.simpleName}.CellPhy.${tree_search_idx}.raxml.bestTree"), path("loglikelihood.${tree_search_idx}.txt")


    script:
    """
    module load cellphy/0.9.2
    raxml-ng-cellphy-linux \
        --search \
        --seed \$RANDOM \
        --msa ${phylo_vcf} \
        --model GTGTR4+G+FO \
        --msa-format VCF \
        --threads ${task.cpus} \
        --prefix ${phylo_vcf.simpleName}.CellPhy.${tree_search_idx} \
        --tree rand{1} \

    loglikelihood=\$(grep "Final LogLikelihood" ${phylo_vcf.simpleName}.CellPhy.${tree_search_idx}.raxml.log | awk '{print \$3}')
    echo \$loglikelihood > loglikelihood.${tree_search_idx}.txt

    """
    stub:
    """
    touch ${phylo_vcf.simpleName}.CellPhy.${tree_search_idx}.raxml.bestTree
    awk -v seed=\$RANDOM 'BEGIN{srand(seed);print -rand()}' > loglikelihood.${tree_search_idx}.txt
    """

}

process BootstrapsCellPhy {
    if ("${workflow.stubRun}" == "false") {
        memory '4 GB'
        cpus 4
    }
    tag "tree-validation"

    publishDir "${params.out}/cellphy/bootstraps", mode: 'symlink'

    input:
    tuple path(phylo_vcf), path(best_tree), val(bootstrap_search_idx)

    output:
    path("${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx}.raxml.bootstraps")


    script:
    """
    module load cellphy/0.9.2
    raxml-ng-cellphy-linux \
        --bootstrap \
        --seed \$RANDOM \
        --msa ${phylo_vcf} \
        --model GTGTR4+G+FO \
        --msa-format VCF \
        --threads ${task.cpus} \
        --prefix ${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx} \
        --bs-trees ${params.bs_trees_per_job} \
        --bs-metric tbe,fbp \


    """
    stub:
    """
    touch ${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx}.raxml.bootstraps
    """

}

process SupportCellPhy {
    if ("${workflow.stubRun}" == "false") {
        memory '8 GB'
        cpus 4
    }
    tag "tree-support"

    publishDir "${params.out}/cellphy/support", mode: 'symlink'

    input:
    path(best_tree)
    path(all_bootstraps)

    output:
    path("${best_tree.simpleName}.CellPhy.raxml.support")


    script:
    """
    module load cellphy/0.9.2

    raxml-ng-cellphy-linux \
        --support \
        --threads ${task.cpus} \
        --tree ${best_tree} \
        --prefix ${best_tree.simpleName}.CellPhy \
        --bs-trees ${all_bootstraps} \


    """
    stub:
    """
    touch ${best_tree.simpleName}.CellPhy.raxml.support
    """

}