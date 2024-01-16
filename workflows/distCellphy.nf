include { CellPhySingleML } from '../modules/phylo'

workflow {
    Channel
        .fromPath( params.joint_vcf )
        .set  { joint_vcf }
    Channel
        .of( 1..params.n_tree_search )
        .set { tree_search_idx }
    CellPhySingleML( joint_vcf, tree_search_idx )
    CellPhySingleML
        .out
        .map { tree, tree_ll -> [ tree, tree_ll.text.toFloat() ]}
        .toSortedList { a, b -> b[1] <=> a[1] }
        .set { best_tree }

    best_tree.view { "The best tree is: ${it[0]}" }
}