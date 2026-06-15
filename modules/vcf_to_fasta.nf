/*
========================================================================================
    MODULE : VCF_TO_FASTA — Conversion VCF multi-sample → alignement FASTA
========================================================================================
    Transforme le VCF multi-sample (bcftools merge) en alignement FASTA
    multi-séquences exploitable par IQ-TREE et PHYLO_REPORT.

    Principe :
      - Chaque souche devient une séquence FASTA
      - Chaque position SNP variable entre les souches devient une colonne
      - Les génotypes manquants (./.) sont convertis en 'N'
      - Les positions non variables (toutes les souches identiques) sont exclues
        car elles n'apportent pas d'information phylogénétique

    Le FASTA produit est fonctionnellement équivalent au snp_alignment.fasta
    produit par kSNP4 en mode de novo — PHYLO_REPORT et IQ-TREE le consomment
    de manière identique dans les deux modes.

    Outil  : vcf_to_fasta.py (script Python maison, bin/)
             Dépendances : stdlib Python uniquement (pas de packages externes)
    Usage  : Mode référence uniquement
    Entrée : merged.vcf.gz (BCFTOOLS_MERGE)
    Sortie : snp_alignment.fasta → IQ-TREE + PHYLO_REPORT
========================================================================================
*/
process VCF_TO_FASTA {
    label 'process_low'
    publishDir "${params.outdir}/snp_alignment", mode: 'copy'

    input:
    // VCF multi-sample fusionné produit par BCFTOOLS_MERGE
    path merged_vcf

    output:
    // Alignement FASTA multi-séquences : une séquence par souche
    path "snp_alignment.fasta", emit: fasta

    script:
    """
    set -euo pipefail

    # projectDir pointe vers la racine du pipeline (où se trouve bin/)
    # Le script utilise uniquement la stdlib Python — compatible avec
    # l'image nicolasdelettre/phyloreport:1.0.1 sans installation supplémentaire
    python3 ${projectDir}/bin/vcf_to_fasta.py \\
        --vcf    ${merged_vcf} \\
        --output snp_alignment.fasta
    """
}
