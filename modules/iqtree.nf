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
    // L'alignement SNP multi-FASTA produit par kSNP4 (snp_alignment.fasta)
    path snp_alignment

    output:
    path "phylo_tree.treefile",                    emit: tree      // Arbre ML principal
    path "phylo_tree.iqtree",                      emit: report    // Rapport détaillé
    path "phylo_tree.log",                         emit: log
    // LEÇON : optional:true → ce fichier n'est pas toujours produit
    // (absent si l'alignement est trop court pour le bootstrap)
    path "phylo_tree.contree", optional: true,     emit: consensus

    script:
    // ─────────────────────────────────────────────────────────────────────────
    // Paramètres IQ-TREE 2 :
    // -s       → fichier d'alignement en entrée
    // -m GTR+G → modèle d'évolution GTR avec variation de taux (bon pour SNPs bactériens)
    // -B 1000  → 1000 réplicats de bootstrap ultra-rapide (UFBoot)
    // -T       → nombre de threads (task.cpus depuis le label)
    // --prefix → préfixe des fichiers de sortie
    // -redo    → force la ré-exécution si un run précédent existe dans le workdir
    //
    // CORRECTION v2.1 : syntaxe modernisée pour IQ-TREE 2.
    //   -nt  → remplacé par -T  (nouveau nom dans IQ-TREE 2)
    //   -bb  → remplacé par -B  (nouveau nom dans IQ-TREE 2)
    // Les anciens noms fonctionnent encore mais sont dépréciés. Autant
    // utiliser la syntaxe moderne puisque notre container est en 2.2.6.
    //
    // LEÇON : Interprétation des valeurs de bootstrap UFBoot :
    //   ≥ 95  → support fort (branche fiable)
    //   70-94 → support modéré
    //   < 70  → support faible (à interpréter avec précaution)
    //
    // LEÇON : pour des alignements SNP (sites variables uniquement), on peut
    // envisager l'ajout de la correction d'ascertainment bias :
    //   -m GTR+ASC+G    → corrige le biais d'échantillonnage des sites invariants
    // À activer si ton arbre montre des branches anormalement longues.
    // ─────────────────────────────────────────────────────────────────────────
    """
    set -euo pipefail

    # Vérification que l'alignement contient assez de séquences
    n_seq=\$(grep -c "^>" ${snp_alignment} || true)
    if [ "\$n_seq" -lt 3 ]; then
        echo "ERREUR : L'alignement contient \$n_seq séquence(s). IQ-TREE requiert au moins 3."
        exit 1
    fi
    echo "Alignement OK : \$n_seq séquences"

    iqtree2 \\
        -s ${snp_alignment} \\
        -m GTR+G \\
        -B 1000 \\
        -T ${task.cpus} \\
        --prefix phylo_tree \\
        -redo
    """
}
