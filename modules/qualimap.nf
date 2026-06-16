/*
========================================================================================
    MODULE : QUALIMAP — Contrôle qualité du mapping BAM
========================================================================================
    Analyse la qualité de l'alignement BAM et génère un rapport HTML par souche.
    Vérifie notamment la couverture du génome, la profondeur de séquençage,
    et la distribution des tailles de reads mappés.

    Deux rapports sont produits :
      - ${sample_id}_qualimap/         : dossier complet lu par MultiQC
                                         (fichiers internes à noms fixes imposés
                                         par Qualimap : qualimapReport.html,
                                         genome_results.txt, etc.)
      - ${sample_id}_qualimapReport.html : copie HTML nommée avec le sample_id
                                           pour identification directe dans EPI2ME
                                           (évite d'avoir tous les rapports
                                           nommés "qualimapReport.html")

    Outil  : Qualimap bamqc
    Usage  : Commun aux deux modes
             Mode de novo    → consomme le BAM de MEDAKA
             Mode référence  → consomme le BAM de MINIMAP2
========================================================================================
*/
process QUALIMAP {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.resultsdir}/qualimap", mode: 'copy'

    input:
    // BAM trié + index : produit par MEDAKA (mode de novo)
    //                    ou MINIMAP2 (mode référence)
    tuple val(sample_id), path(bam), path(bai)

    output:
    // Dossier complet pour MultiQC (noms de fichiers internes imposés par Qualimap)
    path "${sample_id}_qualimap/",           emit: results
    // Copie HTML nommée avec le sample_id pour EPI2ME
    path "${sample_id}_qualimapReport.html", emit: html

    script:
    """
    set -euo pipefail

    qualimap bamqc \\
        -bam ${bam} \\
        -outdir ${sample_id}_qualimap \\
        -outformat HTML \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}

    # Copie du rapport HTML avec le sample_id dans le nom
    # Le dossier ${sample_id}_qualimap/ reste intact pour MultiQC
    cp ${sample_id}_qualimap/qualimapReport.html ${sample_id}_qualimapReport.html
    """
}
