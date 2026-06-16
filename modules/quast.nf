/*
========================================================================================
    MODULE : QUAST — Évaluation de la qualité des assemblages
========================================================================================
    Évalue la qualité structurelle de chaque assemblage de novo :
      - N50, L50 : métriques de contiguïté
      - Nombre de contigs et taille totale
      - Taille du plus grand contig
      - % génome assemblé

    Complémentaire à Qualimap (qualité du mapping) et CheckM2 (complétude biologique).
    Nativement supporté par MultiQC.

    Outil  : QUAST
    Entrée : assemblage polishy par Medaka
    Sortie : rapport TSV/HTML par souche → MultiQC
========================================================================================
*/
process QUAST {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.resultsdir}/quast/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}_quast/", emit: results

    script:
    """
    set -euo pipefail

    quast.py \\
        --output-dir ${sample_id}_quast \\
        --threads ${task.cpus} \\
        --label ${sample_id} \\
        ${assembly}
    """
}
