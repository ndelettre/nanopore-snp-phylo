/*
========================================================================================
    MODULE : BRACKEN — Réestimation des abondances au rang espèce
========================================================================================
    Kraken2 assigne chaque read au plus petit ancêtre commun (LCA) de ses k-mers :
    les reads partagés entre espèces proches restent au rang genre, ce qui fait
    chuter le pourcentage de l'espèce dominante (ex. 45%).
    Bracken redistribue ces reads vers les espèces, de façon probabiliste, à
    partir des distributions de k-mers précalculées de la base Kraken2.

    Outil  : Bracken
    Entrée : rapport Kraken2 (KRAKEN2) + base Kraken2 contenant les fichiers
             databaseXXXmers.kmer_distrib (fournis avec les bases préconstruites)
    Sortie : abondances au rang espèce (TSV) + rapport au format Kraken2
             → le TSV est agrégé dans PHYLO_REPORT
========================================================================================
*/
process BRACKEN {
    tag "${sample_id}"
    publishDir "${params.resultsdir}/bracken", mode: 'copy'

    input:
    // Rapport Kraken2 de la souche
    tuple val(sample_id), path(kraken_report)
    // Base de données Kraken2 (doit contenir les fichiers .kmer_distrib)
    path db

    output:
    // TSV : name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads
    //       new_est_reads  fraction_total_reads
    tuple val(sample_id), path("${sample_id}.bracken.tsv"),    emit: abundance
    // Rapport au format Kraken2 après redistribution (traçabilité)
    tuple val(sample_id), path("${sample_id}.bracken.report"), emit: report

    script:
    """
    set -euo pipefail

    # Les reads Nanopore font plusieurs kb : on prend la plus grande longueur
    # de read pour laquelle la base fournit une distribution de k-mers
    # (300 pour les bases préconstruites PlusPF).
    read_len=\$(ls ${db}/ | sed -n 's/^database\\([0-9]*\\)mers\\.kmer_distrib\$/\\1/p' | sort -n | tail -1)
    if [ -z "\$read_len" ]; then
        echo "ERREUR : aucun fichier databaseXXXmers.kmer_distrib dans ${db}" >&2
        echo "Bracken a besoin des distributions de k-mers (bracken-build)." >&2
        exit 1
    fi
    echo "Bracken : distribution de k-mers pour des reads de \${read_len} pb"

    # -l S  → réestimation au rang espèce
    # -t 10 → seuil minimal de reads par taxon (défaut Bracken)
    bracken \\
        -d ${db} \\
        -i ${kraken_report} \\
        -o ${sample_id}.bracken.tsv \\
        -w ${sample_id}.bracken.report \\
        -r \$read_len \\
        -l S \\
        -t 10
    """
}
