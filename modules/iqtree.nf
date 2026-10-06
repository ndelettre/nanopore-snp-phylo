/*
========================================================================================
    MODULE : IQ-TREE — Construction d'arbre phylogénétique Maximum Likelihood
========================================================================================
    IQ-TREE est l'outil de référence pour les arbres ML en phylogénomique.
    Il est lancé une fois par alignement SNP produit par kSNP4 (tous les SNPs,
    SNPs core) : chaque rapport a ainsi son propre arbre, cohérent avec sa
    matrice de distances.
========================================================================================
*/

process IQTREE {
    tag "${report_name}"
    label 'process_high'
    publishDir "${params.resultsdir}/iqtree", mode: 'copy'

    input:
    // report_name : report_all_snps ou report_core_snps
    tuple val(report_name), path(alignment)

    output:
    tuple val(report_name), path("${report_name}.treefile"), emit: tree
    path "${report_name}.*",                                 emit: all

    script:
    """
    set -euo pipefail

    n_seq=\$(grep -c "^>" ${alignment} || true)
    if [ "\$n_seq" -lt 3 ]; then
        echo "ERREUR : L'alignement contient \$n_seq séquence(s). IQ-TREE requiert au moins 3."
        exit 1
    fi
    echo "Alignement OK : \$n_seq séquences"

    # Ultrafast bootstrap (-B 1000) : IQ-TREE le refuse en dessous de
    # 4 séquences → on construit alors l'arbre sans bootstrap.
    bootstrap_opt=""
    if [ "${params.bootstrap}" = "true" ]; then
        if [ "\$n_seq" -ge 4 ]; then
            bootstrap_opt="-B 1000"
        else
            echo "ATTENTION : moins de 4 séquences, bootstrap désactivé"
        fi
    fi

    # -m GTR+G+ASC : ASC = correction du biais d'échantillonnage, obligatoire
    # pour un alignement ne contenant que des sites variables (SNPs).
    iqtree2 \\
        -s ${alignment} \\
        -m GTR+G+ASC \\
        \$bootstrap_opt \\
        -T ${task.cpus} \\
        --prefix ${report_name} \\
        -redo
    """
}
