/*
========================================================================================
    MODULE : BCFTOOLS_MERGE — Fusion des VCF individuels en VCF multi-sample
========================================================================================
    Combine les VCF produits par Clair3 (un par souche) en un seul fichier
    VCF multi-sample aligné sur les mêmes positions génomiques.

    Les positions présentes dans certaines souches mais absentes d'autres
    sont remplies avec le génotype manquant (./.) — ce comportement est
    géré en aval par vcf_to_fasta.py qui les convertit en 'N'.

    Outil  : bcftools merge
    Usage  : Mode référence uniquement — reçoit TOUS les VCF collectés
    Entrée : Liste de VCF.gz individuels (un par souche)
    Sortie : VCF.gz multi-sample → consommé par VCF_TO_FASTA
========================================================================================
*/
process BCFTOOLS_MERGE {
    label 'process_low'
    publishDir "${params.outdir}/bcftools", mode: 'copy'

    input:
    // Liste de tous les VCF.gz produits par CLAIR3 (.collect() dans main.nf)
    path vcf_list

    output:
    path "merged.vcf.gz", emit: merged_vcf

    script:
    """
    set -euo pipefail

    # Indexation TBI (tabix) de chaque VCF : requis par bcftools merge
    # pour l'accès aléatoire aux positions génomiques
    # Changement de du header de chaque fichier vcf
    for vcf in ${vcf_list}; do
        sample=\$(basename \$vcf .vcf.gz)
        echo "\$sample" > \${sample}_name.txt
        bcftools reheader --samples \${sample}_name.txt \$vcf -o \${sample}_renamed.vcf.gz
        bcftools index --tbi \${sample}_renamed.vcf.gz
    done

    # Fusion multi-sample
    # --merge snps    : stratégie de fusion pour les SNPs
    #                   (résout les conflits aux positions multi-alléliques)
    # --output-type z : sortie compressée BGZF (.vcf.gz)
    bcftools merge \\
        --merge snps \\
        --output-type z \\
        --output merged.vcf.gz \\
        *_renamed.vcf.gz
    """
}
