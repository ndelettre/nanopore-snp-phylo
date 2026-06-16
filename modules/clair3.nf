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
    path ref_fai

    output:
    // VCF.gz par souche, nommé avec le sample_id pour BCFTOOLS_MERGE
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail

    # Variant calling Clair3
    # --platform=ont          : profil d'erreur Oxford Nanopore
    # --model_path            : modèle deep learning (défini dans params.clair3_model)
    # --include_all_ctgs      : appeler les variants sur tous les contigs
    #                           (pas uniquement chr1-22/X/Y comme en mode humain)
    # --no_phasing_for_fa     : désactive le phasage haplotype, non pertinent
    #                           pour les bactéries haploïdes
    # Clair3 résout le .fai depuis le chemin absolu de la référence
    # On crée un lien symbolique avec le bon nom dans le workdir
    REF=\$(realpath ${reference})
    FAI=\$(realpath ${ref_fai})

    # S'assurer que le .fai est à côté de la référence dans le workdir
    ln -sf \${FAI} \${REF}.fai 2>/dev/null || true

run_clair3.sh \\
    --bam_fn=\$(realpath ${bam}) \\
    --ref_fn=\${REF} \\
    --threads=${task.cpus} \\
    --platform=ont \\
    --model_path=/opt/models/${params.clair3_model} \\
    --output=clair3_out \\
    --include_all_ctgs \\
    --no_phasing_for_fa
    --haploid_precise

cp clair3_out/merge_output.vcf.gz ${sample_id}.vcf.gz
    """
}
