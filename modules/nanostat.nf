/*
========================================================================================
    MODULE : NANOSTAT — Statistiques qualité des reads filtrés
========================================================================================
    NanoStat génère un rapport de qualité détaillé sur les reads Nanopore
    (longueur, qualité, N50, total bases...).

    NOTE v2.2 : MultiQC n'a PAS de module natif pour NanoStat. Les fichiers
    seront affichés via le mécanisme "custom content" de MultiQC mais sans
    graphes interactifs dédiés. Pour des graphes natifs, envisager NanoPlot
    ou pycoQC comme alternatives.

    Sorties :
    - {sample_id}_nanostat.txt → stats complètes (longueur, qualité, N50...)
========================================================================================
*/

process NANOSTAT {

    tag "${sample_id}"
    label 'process_mono'

    publishDir "${params.outdir}/nanostat/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    // LEÇON : Pas de tuple ici car MultiQC n'a besoin que du fichier,
    // pas du sample_id. On émet directement le path.
    path "${sample_id}_nanostat.txt", emit: stats

    script:
    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : task.cpus est une variable Nextflow automatique.
    // Elle prend la valeur du label associé au process dans nextflow.config.
    // Ici process_low → cpus = 2. C'est plus fiable qu'un params.threads
    // global qui peut diverger des ressources réellement allouées.
    // ─────────────────────────────────────────────────────────────────────────
    """
    set -euo pipefail

    NanoStat \\
        --fastq ${fastq} \\
        --name "${sample_id}_nanostat.txt" \\
        --threads ${task.cpus}
    """
}
