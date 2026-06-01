/*
========================================================================================
    MODULE : IQ-TREE — Construction d'arbre phylogénétique Maximum Likelihood
========================================================================================
    IQ-TREE est l'outil de référence pour les arbres ML en phylogénomique.
    Il prend en entrée l'alignement SNP FASTA produit par kSNP4
    et produit un arbre au format Newick (.treefile).
========================================================================================
*/

process IQTREE {
    label 'process_high'
    publishDir "${params.outdir}/iqtree", mode: 'copy'

    input:
    path alignment

    output:
    path "phylo_tree.treefile", emit: tree
    path "phylo_tree.*",        emit: all

    script:
    def bootstrap = params.bootstrap ? "-B 1000" : ""
    """
    set -euo pipefail

    n_seq=\$(grep -c "^>" snp_alignment.fasta || true)
    if [ "\$n_seq" -lt 3 ]; then
        echo "ERREUR : L'alignement contient \$n_seq séquence(s). IQ-TREE requiert au moins 3."
        exit 1
    fi
    echo "Alignement OK : \$n_seq séquences"

    iqtree2 \\
        -s snp_alignment.fasta \\
        -m GTR+G \\
        ${bootstrap} \\
        -T 12 \\
        --prefix phylo_tree \\
        -redo
    """
}
