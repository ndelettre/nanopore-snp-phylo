/*
========================================================================================
    MODULE : CLAIR3 — Variant calling par deep learning sur reads longs ONT
========================================================================================
    Détecte les SNPs et indels par souche en alignant les reads sur la référence.
    Utilise une approche en deux passes :
      1. Pileup model   : rapide, détecte les candidats variants
      2. Full-alignment : précis, confirme et filtre les candidats

    Outil  : Clair3 (HKU-BAL)
    Usage  : Mode référence uniquement — une instance par souche
    Entrée : BAM trié+indexé (Minimap2) + référence FASTA
    Sortie : VCF.gz par souche → fusionné ensuite par BCFTOOLS_MERGE
========================================================================================
*/
process CLAIR3 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/clair3/${sample_id}", mode: 'copy'

    input:
    // BAM trié + index produits par MINIMAP2
    tuple val(sample_id), path(bam), path(bai)
    // Génome de référence FASTA (partagé entre toutes les souches)
    path reference

    output:
    // VCF.gz par souche, nommé avec le sample_id pour BCFTOOLS_MERGE
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail

    # L'index FASTA est requis par Clair3 pour l'accès aléatoire aux contigs
    samtools faidx ${reference}

    # Variant calling Clair3
    # --platform=ont          : profil d'erreur Oxford Nanopore
    # --model_path            : modèle deep learning (défini dans params.clair3_model)
    # --include_all_ctgs      : appeler les variants sur tous les contigs
    #                           (pas uniquement chr1-22/X/Y comme en mode humain)
    # --no_phasing_for_fa     : désactive le phasage haplotype, non pertinent
    #                           pour les bactéries haploïdes
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=/opt/models/${params.clair3_model} \\
        --output=clair3_out \\
        --include_all_ctgs \\
        --no_phasing_for_fa

    # Renommer la sortie avec le sample_id pour traçabilité
    # merge_output.vcf.gz = VCF final consolidé produit par Clair3
    cp clair3_out/merge_output.vcf.gz ${sample_id}.vcf.gz
    """
}
