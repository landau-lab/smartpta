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
        --msa ${phylo_vcf} \
        --model GTGTR4+G+FO \
        --msa-format VCF \
        --threads ${task.cpus} \
        --prefix ${phylo_vcf.simpleName}.CellPhy.${tree_search_idx} \
        --tree rand{1} \

    loglikelihood=\$(grep "Final LogLikelihood" ${phylo_vcf.simpleName}.CellPhy.raxml.log | awk '{print \$3}')
    echo \$loglikelihood > loglikelihood.${tree_search_idx}.txt

    """
    stub:
    """
    touch ${phylo_vcf.simpleName}.CellPhy.${tree_search_idx}.raxml.bestTree
    awk -v seed=\$RANDOM 'BEGIN{srand(seed);print -rand()}' > loglikelihood.${tree_search_idx}.txt
    """

}

process CellPhyBootstraps {
    if ("${workflow.stubRun}" == "false") {
        memory '512 GB'
        cpus 36
        queue 'bigmem'
    }
    tag "tree-validation"

    publishDir "${params.out}/cellphy/bootstraps", mode: 'symlink'

    input:
    path(phylo_vcf)
    path(best_tree)
    each bootstrap_search_idx

    output:
    path("${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx}.raxml.support")


    script:
    """
    module load cellphy/0.9.2
    raxml-ng-cellphy-linux \
        --bootstrap \
        --msa ${phylo_vcf} \
        --model GTGTR4+G+FO \
        --msa-format VCF \
        --threads ${task.cpus} \
        --prefix ${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx} \
        --bs-trees ${params.bs_trees_per_job} \
        --bs-metric tbe,fbp \

    raxml-ng-cellphy-linux \
        --support \
        --threads ${task.cpus} \
        --tree ${best_tree} \
        --bs-trees ${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx}.raxml.bootstraps \


    """
    stub:
    """
    touch ${phylo_vcf.simpleName}.CellPhy.${bootstrap_search_idx}.raxml.support
    """

}