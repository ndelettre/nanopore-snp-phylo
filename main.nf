#!/usr/bin/env nextflow
/*
========================================================================================
    PIPELINE NANOPORE - ANALYSE DE SOUCHES BACTÉRIENNES v2.4
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    WORKFLOW :
       NANOFILT → NANOSTAT → FLYE → MEDAKA → QUALIMAP
                                           → KSNP4 → IQTREE → PHYLO_REPORT
                                           → MULTIQC
========================================================================================
*/

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMÈTRES DU PIPELINE
// ─────────────────────────────────────────────────────────────────────────────
params.fastq_dir    = null          // Dossier contenant les FASTQ (obligatoire)
params.outdir       = "output"      // Dossier de sortie des rapports
params.resultsdir   = "results"     // Dossier de sortie des fichiers
params.min_length   = 200           // Longueur minimale des reads (NanoFilt)
params.min_quality  = 10            // Qualité minimale des reads Q-score (NanoFilt)
params.genome_size  = "5m"          // Taille estimée du génome pour Flye (ex: 5m = 5 Mb)
params.medaka_model = "r1041_e82_400bps_sup_v5.2.0"
                                    // Modèle Medaka : r1041 = R10.4.1 | e82 = Kit 14 | sup = SUP
params.bootstrap    = true          // Active le calcul des valeurs de bootstrap IQ-TREE
params.kraken_db = "/data/kraken2_db"

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// ─────────────────────────────────────────────────────────────────────────────
include { NANOFILT }     from './modules/nanofilt.nf'
include { NANOSTAT }     from './modules/nanostat.nf'
include { FLYE }         from './modules/flye.nf'
include { MEDAKA }       from './modules/medaka.nf'
include { QUALIMAP }     from './modules/qualimap.nf'
include { KSNP4 }        from './modules/ksnp4.nf'
include { IQTREE }       from './modules/iqtree.nf'
include { MULTIQC }      from './modules/multiqc.nf'
include { PHYLO_REPORT } from './modules/phylo_report.nf'
include { KRAKEN2 } from './modules/kraken2.nf'

// ─────────────────────────────────────────────────────────────────────────────
// BANNIÈRE DE DÉMARRAGE
// ─────────────────────────────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════╗
║     PIPELINE SNP & PHYLOGÉNIE - NANOPORE MINION  v2.4    ║
╚══════════════════════════════════════════════════════════╝
  Dossier FASTQ   : ${params.fastq_dir}
  Dossier sortie  : ${params.outdir}
  Qualité min.    : ${params.min_quality}
  Longueur min.   : ${params.min_length} bp
  Taille génome   : ${params.genome_size}
  Modèle Medaka   : ${params.medaka_model}
  Bootstrap       : ${params.bootstrap}
──────────────────────────────────────────────────────────
""".stripIndent()

// ─────────────────────────────────────────────────────────────────────────────
// VALIDATION DES PARAMÈTRES OBLIGATOIRES
// ─────────────────────────────────────────────────────────────────────────────
if (!params.fastq_dir) {
    error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
}
if (!params.kraken_db) {
    error "ERREUR : --kraken_db est obligatoire.\nUsage : nextflow run main.nf --kraken_db /chemin/vers/db"
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Création du channel d'entrée ─────────────────────────────────────────
    ch_fastq = Channel
        .fromPath(
            "${params.fastq_dir}/*.{fastq,fastq.gz,fq,fq.gz}",
            checkIfExists: true
        )
        .map { file ->
            def sample_id = file.getSimpleName()
                .replaceAll(/\.fastq.*/, '')
                .replaceAll(/\.fq.*/, '')
            tuple(sample_id, file)
        }

    // ── Filtrage et QC des reads ──────────────────────────────────────────────
    NANOFILT(ch_fastq)
    NANOSTAT(NANOFILT.out.reads)

    //vérification de l'identification
    ch_kraken_db = Channel.fromPath(params.kraken_db, checkIfExists: true).first()
    KRAKEN2(NANOFILT.out.reads, ch_kraken_db)

    // ── Assemblage de novo ────────────────────────────────────────────────────
    // Flye est optimisé pour les reads longs avec taux d'erreur élevé
    FLYE(NANOFILT.out.reads)

    // ── Polissage ─────────────────────────────────────────────────────────────
    // .join() garantit que chaque souche reçoit SES propres reads et assemblage
    ch_medaka_input = NANOFILT.out.reads.join(FLYE.out.assembly)
    MEDAKA(ch_medaka_input)

    // ── QC du mapping ─────────────────────────────────────────────────────────
    QUALIMAP(MEDAKA.out.bam)

    // ── SNP calling multi-souches ─────────────────────────────────────────────
    // kSNP4 compare tous les assemblages simultanément via approche k-mer
    // .collect() attend que toutes les souches soient assemblées
    ch_all_assemblies = MEDAKA.out.assembly
        .map { sample_id, fasta -> fasta }
        .collect()

    KSNP4(ch_all_assemblies)

    // ── Arbre phylogénétique ──────────────────────────────────────────────────
    IQTREE(KSNP4.out.snp_alignment)

    // ── Rapports phylogénétiques interactifs ──────────────────────────────────
    // Deux rapports : tous les SNPs et SNPs core uniquement
    ch_all_snps = KSNP4.out.snp_alignment
        .map { fasta -> tuple("report_all_snps", fasta) }
        .combine(IQTREE.out.tree)

    ch_core_snps = KSNP4.out.core_snp_matrix
        .map { fasta -> tuple("report_core_snps", fasta) }
        .combine(IQTREE.out.tree)

    PHYLO_REPORT(ch_all_snps.mix(ch_core_snps))

    // ── Rapport MultiQC ───────────────────────────────────────────────────────
    ch_multiqc_files = NANOSTAT.out.stats
        .mix(FLYE.out.stats)
        .mix(QUALIMAP.out.results)
        .collect()

    MULTIQC(ch_multiqc_files)
}

// ─────────────────────────────────────────────────────────────────────────────
// RÉSUMÉ DE FIN DE PIPELINE
// ─────────────────────────────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
╔══════════════════════════════════════════════════════════╗
║                  PIPELINE TERMINÉ !                      ║
╚══════════════════════════════════════════════════════════╝
  Statut    : ${workflow.success ? '✅ Succès' : '❌ Échec'}
  Durée     : ${workflow.duration}
  Rapports : ${params.outdir}/
  Résultats : ${params.resultsdir}/
──────────────────────────────────────────────────────────
""".stripIndent()
}
