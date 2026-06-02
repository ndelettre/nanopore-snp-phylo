/*
========================================================================================
    MODULE : SNIPPY — SNP calling avec génome de référence (mode optionnel)
========================================================================================
    Snippy identifie les variants (SNPs, indels) de chaque souche
    par rapport à un génome de référence connu.
    Il est activé uniquement si --reference est fourni en ligne de commande.
========================================================================================
*/

process SNIPPY {

    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/snippy/${sample_id}", mode: 'copy'

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : .combine() dans main.nf a produit le tuple :
    //   [sample_id, assembly_polished.fasta, reference.fasta]
    //
    // On utilise les assemblages polishés Medaka (FASTA) et non les reads
    // bruts FASTQ, car l'option --contigs de Snippy attend des contigs assemblés.
    // Pour des reads non assemblés, il faudrait utiliser --se (single-end).
    // ─────────────────────────────────────────────────────────────────────────
    input:
    tuple val(sample_id), path(polished_assembly), path(reference)

    // ─────────────────────────────────────────────────────────────────────────
    // CORRECTION v2.2 : `optional: true` ne peut PAS être combiné avec
    // `tuple val(...) path(...)` de la manière écrite en v2.1 — Nextflow
    // rejette cette syntaxe. On sépare l'émission optionnelle en `path` pur
    // (sans sample_id dans le tuple). Si besoin du sample_id en aval, on le
    // récupère via un .map() dans le workflow.
    //
    // Snippy produit par défaut {prefix}.vcf dans le dossier --outdir.
    // Ici prefix = ${sample_id}, outdir = ${sample_id}_snippy.
    // Donc le VCF final est : ${sample_id}_snippy/${sample_id}.vcf
    // ─────────────────────────────────────────────────────────────────────────
    output:
    tuple val(sample_id), path("${sample_id}_snippy/"),                 emit: results
    tuple val(sample_id), path("${sample_id}_snippy/${sample_id}.vcf"), emit: vcf
    path "${sample_id}_snippy/${sample_id}.aligned.fa", optional: true, emit: aligned

    script:
    // ─────────────────────────────────────────────────────────────────────────
    // CORRECTION v2.1 : `--ctgs` est déprécié depuis Snippy 4.4.
    // Nouvelle option : `--contigs`. Fonctionnellement identique.
    //
    // Paramètres Snippy :
    // --ref      → génome de référence au format FASTA ou GenBank
    // --contigs  → contigs de l'assemblage polished (FASTA)
    // --outdir   → dossier de sortie
    // --prefix   → préfixe des fichiers de sortie
    // --cpus     → nombre de threads
    // --force    → écrase le dossier de sortie s'il existe (utile pour -resume)
    // ─────────────────────────────────────────────────────────────────────────
    """
    set -euo pipefail

    snippy \\
        --cpus ${task.cpus} \\
        --outdir ${sample_id}_snippy \\
        --ref ${reference} \\
        --contigs ${polished_assembly} \\
        --prefix ${sample_id} \\
        --force

    # Vérification que le VCF a bien été produit
    if [ ! -f "${sample_id}_snippy/${sample_id}.vcf" ]; then
        echo "ERREUR : Snippy n'a pas produit le VCF attendu pour ${sample_id}"
        ls -la ${sample_id}_snippy/ || true
        exit 1
    fi
    """
}
