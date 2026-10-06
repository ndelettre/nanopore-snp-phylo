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
    Entrée : reads filtrés (CHOPPER) + base de données PlusPF-8
    Sortie : rapport taxonomique par souche (format Kraken2 standard)
             → agrégé dans PHYLO_REPORT pour le tableau HTML
========================================================================================
*/
process KRAKEN2 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.resultsdir}/kraken2", mode: 'copy'

    input:
    // Reads filtrés par Chopper
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

/*
    KRAKEN2_REFERENCE — Espèce de la souche de référence (--reference_fasta)
    La référence est un assemblage : on classe ses contigs et on garde la
    sortie par séquence (taxid + longueur). PHYLO_REPORT retient l'espèce qui
    couvre le plus de paires de bases, pour qu'un plasmide classé ailleurs ne
    l'emporte pas sur le chromosome. Pas de Bracken ici : il travaille sur des
    reads, pas sur quelques contigs.
*/
process KRAKEN2_REFERENCE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.resultsdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path db

    output:
    tuple val(sample_id),
          path("${sample_id}.ref_kraken2.report"),
          path("${sample_id}.ref_kraken2.out"),     emit: classification

    script:
    """
    set -euo pipefail

    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --report ${sample_id}.ref_kraken2.report \\
        --output ${sample_id}.ref_kraken2.out \\
        ${assembly}
    """
}
