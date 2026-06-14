/*
========================================================================================
    MODULE : MULTIQC — Rapport de qualité agrégé
========================================================================================
    MultiQC agrège les rapports de qualité de tous les outils en un seul
    rapport HTML interactif. Il détecte automatiquement les formats compatibles
    dans le dossier courant.

    Sources parsées ici :
    - NanoStat (.txt) → stats qualité des reads filtrés par échantillon
    - Flye (assembly_info.txt) → stats des assemblages (N50, taille totale...)
========================================================================================
*/

process MULTIQC {

    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: 'copy'
    publishDir "${params.outdir}/output",         mode: 'copy', pattern: "multiqc_report.html"

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : .mix() dans main.nf a fusionné les channels NanoStat et Flye.
    // .collect() a ensuite attendu que tous les fichiers soient disponibles.
    // MultiQC reçoit donc une liste plate de tous les fichiers de stats.
    // ─────────────────────────────────────────────────────────────────────────
    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data/",       emit: data

    script:
    """
    set -euo pipefail

    multiqc . \\
        --filename multiqc_report.html \\
        --title "Rapport QC - Pipeline Nanopore SNP" \\
        --force
    """
}
