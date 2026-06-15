/*
========================================================================================
    MODULE : MINIMAP2 — Mapping des reads longs sur une référence
========================================================================================
    Aligne les reads ONT filtrés sur le génome de référence fourni.
    Produit un BAM trié et indexé par échantillon.

    Outil    : Minimap2 + Samtools (inclus dans le container)
    Preset   : map-ont → optimisé pour les reads longs Oxford Nanopore
               (gestion des taux d'erreur élevés et des indels fréquents)
    Usage    : Mode référence uniquement
========================================================================================
*/
process MINIMAP2 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/minimap2", mode: 'copy'

    input:
    // Tuple : identifiant de l'échantillon + fichier FASTQ filtré par NanoFilt
    tuple val(sample_id), path(reads)
    // Génome de référence FASTA (partagé entre toutes les souches)
    path reference

    output:
    // BAM trié + index .bai : nécessaires pour Clair3 et Qualimap
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    set -euo pipefail

    # Alignement des reads sur la référence, tri et indexation en une seule commande
    # -ax map-ont : preset Oxford Nanopore long reads
    # -t          : nombre de threads alloués par Nextflow selon le label
    # Le pipe évite un fichier SAM intermédiaire (gain d'espace disque)
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}
