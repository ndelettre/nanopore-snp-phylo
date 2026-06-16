/*
========================================================================================
    MODULE : CHECKM2 — Évaluation de la qualité des assemblages bactériens
========================================================================================
    Évalue la complétude et la contamination de chaque assemblage via un modèle
    de deep learning entraîné sur des génomes bactériens de référence.

    Seuils qualité standard (MIMAG) :
      - Haute qualité  : complétude ≥90% ET contamination <5%
      - Moyenne qualité : complétude ≥50% ET contamination <10%
      - En dessous     : assemblage à exclure de l'analyse phylogénétique

    Outil  : CheckM2
    Entrée : assemblages polishés (MEDAKA) — un par souche
    Sortie : rapport TSV avec complétude + contamination par souche
             → intégré dans PHYLO_REPORT avec alertes visuelles
========================================================================================
*/
process CHECKM2 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.resultsdir}/checkm2", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_checkm2/quality_report.tsv"), emit: report

    script:
    """
    set -euo pipefail

    checkm2 predict \\
        --input ${assembly} \\
        --output-directory ${sample_id}_checkm2 \\
        --threads ${task.cpus} \\
        --force
    """
}
