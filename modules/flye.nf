/*
========================================================================================
    MODULE : FLYE — Assemblage de novo pour lectures longues Nanopore
========================================================================================
*/

process FLYE {

    tag "${sample_id}"
    label 'process_medium'   // Flye est gourmand en CPU/RAM → label 'high'

    publishDir "${params.outdir}/flye/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_fastq)

    output:
    // assembly.fasta est le fichier de sortie standard de Flye
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    // Stats d'assemblage (N50, taille totale...) parsées par MultiQC
    path "assembly_info_${sample_id}.txt",                     emit: stats
    // Log complet de Flye pour le débogage
    path "flye_${sample_id}.log",                              emit: log

    script:
    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : "set -euo pipefail" est essentiel ici.
    // Sans "pipefail", si Flye échoue dans un pipe bash, le code de retour
    // serait celui de la dernière commande (souvent 0), masquant l'erreur.
    // Nextflow penserait que tout va bien et continuerait le pipeline.
    //
    // LEÇON : On sépare l'exécution de Flye de la redirection du log
    // ("> log 2>&1" au lieu de "| tee log") pour la même raison.
    //
    // LEÇON : task.cpus reflète automatiquement le label 'process_high'
    // défini dans nextflow.config (cpus = 12 dans notre config).
    //
    // LEÇON : --nano-hq est pour les données récentes (R10.x, Dorado sup).
    //         Utilise --nano-raw pour les anciennes données R9.x.
    // ─────────────────────────────────────────────────────────────────────────
    """
    set -euo pipefail

    flye \\
        --nano-hq ${filtered_fastq} \\
        --out-dir flye_output \\
        --genome-size ${params.genome_size} \\
        --threads ${task.cpus} \\
        > flye_${sample_id}.log 2>&1

    # Vérification que l'assemblage a bien été produit
    if [ ! -f flye_output/assembly.fasta ]; then
        echo "ERREUR : Flye n'a pas produit d'assemblage pour ${sample_id}"
        echo "Contenu du dossier de sortie :"
        ls -la flye_output/ || true
        exit 1
    fi

    # Renommage des fichiers de sortie avec le nom de l'échantillon
    cp flye_output/assembly.fasta ${sample_id}_assembly.fasta
    cp flye_output/assembly_info.txt assembly_info_${sample_id}.txt
    """
}
