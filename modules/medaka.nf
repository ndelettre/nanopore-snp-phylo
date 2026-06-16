/*
========================================================================================
    MODULE : MEDAKA — Polissage d'assemblage Nanopore
========================================================================================
    Medaka corrige les erreurs d'assemblage en remappant les reads originaux
    sur l'assemblage Flye. Il utilise un modèle de réseau de neurones
    entraîné spécifiquement pour chaque chimie de flowcell Nanopore.
========================================================================================
*/

process MEDAKA {

    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.resultsdir}/medaka/${sample_id}", mode: 'copy'

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : Ici l'input est un tuple à 3 éléments.
    // Il provient du .join() dans main.nf qui combine :
    //   - NANOFILT.out.reads    → [sample_id, fastq_filtré]
    //   - FLYE.out.assembly     → [sample_id, assembly.fasta]
    // Le .join() les fusionne en : [sample_id, fastq_filtré, assembly.fasta]
    // ─────────────────────────────────────────────────────────────────────────
    input:
    tuple val(sample_id), path(filtered_fastq), path(draft_assembly)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), emit: assembly
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    // ─────────────────────────────────────────────────────────────────────────
    // NOTE : Choix du modèle Medaka selon ta flowcell :
    //   r941_min_high_g360          → MinION R9.4.1, Guppy high accuracy
    //   r941_min_sup_g507           → MinION R9.4.1, Guppy super accuracy
    //   r1041_e82_400bps_sup_v5.0.0 → MinION R10.4.1, Dorado sup
    // Lance "medaka tools list_models" pour voir tous les modèles disponibles.
    // ─────────────────────────────────────────────────────────────────────────
    """
    set -euo pipefail

    medaka_consensus \\
        -i ${filtered_fastq} \\
        -d ${draft_assembly} \\
        -o medaka_output \\
        -m ${params.medaka_model} \\
        -t ${task.cpus}

    # Vérification que le consensus a bien été produit
    if [ ! -f medaka_output/consensus.fasta ]; then
        echo "ERREUR : Medaka n'a pas produit de consensus pour ${sample_id}"
        exit 1
    fi

    # Renommage du consensus final
    cp medaka_output/consensus.fasta ${sample_id}_polished.fasta

    # Exposer le BAM du polishing pour Qualimap
    cp medaka_output/calls_to_draft.bam     ${sample_id}.bam
    cp medaka_output/calls_to_draft.bam.bai ${sample_id}.bam.bai
    """
}
