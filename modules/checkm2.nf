/*
========================================================================================
    MODULE : CHECKM2 — Évaluation de la qualité des assemblages bactériens
========================================================================================
    Évalue la complétude et la contamination de chaque assemblage via un modèle
    de deep learning entraîné sur des génomes bactériens de référence.

    Seuils appliqués dans le rapport (PHYLO_REPORT) :
      - Complétude    ≥ 99%
      - Contamination < 1%
    Plus stricts que MIMAG (≥90% / <5%), adaptés à des isolats purs
    séquencés en génome complet. Une souche hors seuils est signalée NOK
    dans le rapport, mais n'est pas exclue de la phylogénie.

    Outil  : CheckM2
    Entrée : assemblages polishés (MEDAKA) — un par souche
    Sortie : rapport TSV avec complétude + contamination par souche
             → intégré dans PHYLO_REPORT (statut OK/NOK)
========================================================================================
*/
process CHECKM2 {
    tag "${sample_id}"
    publishDir "${params.resultsdir}/checkm2", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path checkm2_db

    output:
tuple val(sample_id), path("${sample_id}_checkm2_quality_report.tsv"), emit: report

script:
"""
set -euo pipefail

checkm2 predict \\
    --input ${assembly} \\
    --output-directory ${sample_id}_checkm2 \\
    --threads ${task.cpus} \\
    --database_path ${checkm2_db} \\
    --force

cp ${sample_id}_checkm2/quality_report.tsv ${sample_id}_checkm2_quality_report.tsv
"""
}
