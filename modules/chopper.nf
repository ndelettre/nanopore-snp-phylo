/*
========================================================================================
    MODULE : CHOPPER — Filtrage qualité et longueur des reads
========================================================================================
    LEÇON : Un module Nextflow contient un seul "process".
    Un process a toujours cette structure :

        process NOM {
            directive       ← comment exécuter (Docker, CPU, RAM, sortie...)
            input:          ← ce qu'il reçoit
            output:         ← ce qu'il produit
            script:         ← la commande shell à exécuter
        }

    Chopper remplace NanoFilt (même auteur, W. De Coster) : même filtrage
    sur la qualité moyenne et la longueur des reads, mais écrit en Rust,
    multithread et maintenu. Il lit le FASTQ sur l'entrée standard.
========================================================================================
*/

process CHOPPER {

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : Les directives définissent le comportement du process.
    //
    // tag       → affiche le nom de l'échantillon dans les logs (très utile)
    // label     → regroupe les process par profil de ressources (défini dans nextflow.config)
    // publishDir → copie les fichiers de sortie dans le dossier résultats
    //              "mode: 'copy'" copie le fichier (vs 'symlink' qui crée un lien)
    // ─────────────────────────────────────────────────────────────────────────
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.resultsdir}/chopper/${sample_id}", mode: 'copy'

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : L'input est un tuple (paire) : [identifiant, fichier]
    // Le mot-clé "val" reçoit une valeur simple (texte, nombre)
    // Le mot-clé "path" reçoit un fichier (Nextflow le gère automatiquement)
    // ─────────────────────────────────────────────────────────────────────────
    input:
    tuple val(sample_id), path(fastq)

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : L'output déclare ce que le process va produire.
    // On utilise des noms (.reads, .log) pour les référencer dans main.nf
    // avec la syntaxe CHOPPER.out.reads
    //
    // emit: → donne un nom à cet output pour y accéder facilement
    // ─────────────────────────────────────────────────────────────────────────
    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: reads
    path "${sample_id}_chopper.log",                               emit: log

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : "set -euo pipefail" est une bonne pratique bash :
    //   -e  → arrête le script à la première erreur
    //   -u  → erreur si une variable non définie est utilisée
    //   -o pipefail → propage les erreurs à travers les pipes (|)
    //
    // Le FASTQ (compressé ou non) est décompressé à la volée vers chopper,
    // sans fichier intermédiaire sur le disque.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    def reader = fastq.name.endsWith('.gz') ? 'zcat' : 'cat'
    """
    set -euo pipefail

    # Filtrage : qualité moyenne ≥ min_quality et longueur ≥ min_length.
    # chopper écrit son bilan (reads gardés / total) sur stderr → log.
    ${reader} "${fastq}" \\
        | chopper \\
            --quality ${params.min_quality} \\
            --minlength ${params.min_length} \\
            --threads ${task.cpus} \\
            2> "${sample_id}_chopper.log" \\
        | gzip > "${sample_id}_filtered.fastq.gz"

    cat "${sample_id}_chopper.log"

    # Vérification que le fichier filtré n'est pas vide (0 read survivant) :
    # le pipeline s'arrête volontairement si une souche n'a aucun read.
    # NOTE : on désactive temporairement pipefail pour la vérification
    # car `zcat | head -n 4` provoque un SIGPIPE sur zcat quand head ferme
    # le pipe après 4 lignes. Avec `set -o pipefail`, SIGPIPE fait échouer
    # la commande entière et déclenche `set -e` → process killé à tort.
    set +o pipefail
    n_lines=\$(zcat "${sample_id}_filtered.fastq.gz" | head -n 4 | wc -l)
    set -o pipefail
    if [ "\$n_lines" -lt 4 ]; then
        echo "ERREUR : aucun read n'a survécu au filtrage pour ${sample_id}" >&2
        echo "  → vérifier --min_quality (${params.min_quality}) et --min_length (${params.min_length})" >&2
        exit 1
    fi
    """
}
