/*
========================================================================================
    MODULE : NANOFILT
========================================================================================
    LEÇON : Un module Nextflow contient un seul "process".
    Un process a toujours cette structure :

        process NOM {
            directive       ← comment exécuter (Docker, CPU, RAM, sortie...)
            input:          ← ce qu'il reçoit
            output:         ← ce qu'il produit
            script:         ← la commande shell à exécuter
        }
========================================================================================
*/

process NANOFILT {

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

    publishDir "${params.outdir}/nanofilt/${sample_id}", mode: 'copy'

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
    // avec la syntaxe NANOFILT.out.reads
    //
    // emit: → donne un nom à cet output pour y accéder facilement
    // ─────────────────────────────────────────────────────────────────────────
    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: reads
    path "${sample_id}_nanofilt.log",                              emit: log

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : Le script est du code shell classique.
    // Les variables Nextflow (sample_id, params.xxx) s'utilisent avec ${}.
    // On peut écrire du bash multi-lignes ici.
    //
    // LEÇON : "set -euo pipefail" est une bonne pratique bash :
    //   -e  → arrête le script à la première erreur
    //   -u  → erreur si une variable non définie est utilisée
    //   -o pipefail → propage les erreurs à travers les pipes (|)
    //
    // CORRECTION v2.1 (bug BLOQUANT) : la version précédente avait une
    // syntaxe bash invalide. Elle tentait de chaîner un bloc if/else/fi
    // avec un pipe `| NanoFilt ...` sur la ligne suivante :
    //
    //     if [[ ... ]]; then gunzip -c ...; else cat ...; fi
    //     | NanoFilt ...    ← INVALIDE : un pipe ne peut pas démarrer une ligne
    //
    // Solution propre : on utilise une substitution de process bash `< <(...)`
    // qui injecte le stdout du bloc if/else directement dans NanoFilt via stdin.
    // Alternative plus lisible : décompresser d'abord dans un fichier
    // intermédiaire, puis appliquer NanoFilt dessus (approche retenue ici
    // car plus simple à déboguer et compatible avec tous les shells POSIX).
    //
    // CORRECTION v2.1 (robustesse) : on vérifie que le fichier de sortie
    // n'est pas vide après filtrage. Sans cette vérification, si NanoFilt
    // écarte 100% des reads, gzip produit un fichier valide mais vide qui
    // ferait planter Flye en aval avec un message obscur.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    """
    set -euo pipefail

    # Étape 1 : décompression conditionnelle vers un fichier temporaire
    if [[ "${fastq}" == *.gz ]]; then
        gunzip -c "${fastq}" > input_raw.fastq
    else
        cp "${fastq}" input_raw.fastq
    fi

    # Étape 2 : filtrage NanoFilt → gzip en sortie
    NanoFilt \\
        --quality ${params.min_quality} \\
        --length ${params.min_length} \\
        --logfile "${sample_id}_nanofilt.log" \\
        < input_raw.fastq \\
        | gzip > "${sample_id}_filtered.fastq.gz"

    # Nettoyage du fichier temporaire
    rm -f input_raw.fastq

    # Garantit l'existence du log même si NanoFilt ne l'a pas créé
    touch "${sample_id}_nanofilt.log"

    # Vérification que le fichier filtré n'est pas vide (0 read survivant)
    # NOTE v2.2 : on désactive temporairement pipefail pour la vérification
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
