/*
========================================================================================
    MODULE : MLST — Typage multilocus des séquences (Multi-Locus Sequence Typing)
========================================================================================
    Détermine le sequence type (ST) de chaque souche à partir de son assemblage.
    L'outil mlst de Torsten Seemann détecte automatiquement l'espèce et applique
    le schéma MLST approprié depuis la base PubMLST.

    Pour Clostridioides difficile : schéma cdifficile (7 loci : adk, atpA, dxr,
    glyA, recA, sodA, tpi) → ST directement comparable aux bases épidémiologiques.

    Outil  : mlst (Seemann T, github.com/tseemann/mlst)
    Entrée : assemblages polishés (MEDAKA)
    Sortie : fichier TSV par souche avec ST et allèles
             → intégré dans PHYLO_REPORT
========================================================================================
*/
process MLST {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.resultsdir}/mlst", mode: 'copy'

    input:
    // Assemblage polishy par Medaka
    tuple val(sample_id), path(assembly)

    output:
    // TSV : sample  schéma  ST  allèle1  allèle2 ...
    tuple val(sample_id), path("${sample_id}.mlst.tsv"), emit: mlst

    script:
    """
    set -euo pipefail

    mlst \\
        --threads ${task.cpus} \\
        --json ${sample_id}.mlst.json \\
        ${assembly} > ${sample_id}.mlst.tsv
    """
}
