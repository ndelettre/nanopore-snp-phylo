/*
========================================================================================
    MODULE : KRAKEN2 — Identification taxonomique et contrôle de pureté
========================================================================================
    Classifie les reads filtrés contre une base de données Kraken2 (PlusPF-8)
    pour identifier l'espèce dominante et détecter les contaminations.

    Utilisé comme contrôle qualité avant l'assemblage :
      - Confirme l'identité de l'espèce séquencée
      - Détecte les contaminations (autre espèce >1% des reads)
      - Alerte si la souche n'est pas pure

    Outil  : Kraken2
    Entrée : reads filtrés (NANOFILT) + base de données PlusPF-8
    Sortie : rapport taxonomique par souche (format Kraken2 standard)
             → agrégé dans PHYLO_REPORT pour le tableau HTML
========================================================================================
*/
process KRAKEN2 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    // Reads filtrés par NanoFilt
    tuple val(sample_id), path(reads)
    // Base de données Kraken2 (PlusPF-8 recommandée)
    path db

    output:
    // Rapport taxonomique : une ligne par taxon avec pourcentages
    // Format : %reads  nb_reads_clade  nb_reads_taxon  rang  taxid  nom
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report

    script:
    """
    set -euo pipefail

    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --report ${sample_id}.kraken2.report \\
        --output /dev/null \\
        ${reads}
    """
}
