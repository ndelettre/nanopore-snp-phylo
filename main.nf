#!/usr/bin/env nextflow
/*
========================================================================================
    PIPELINE NANOPORE - ANALYSE DE SOUCHES BACTÉRIENNES v2.4
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    WORKFLOW :
       NANOFILT → NANOSTAT → KRAKEN2 (identification espèce)
                          → FLYE → MEDAKA → QUALIMAP
                                          → QUAST
                                          → MLST
                                          → CHECKM2
                                          → KSNP4 → IQTREE → PHYLO_REPORT
                                          → MULTIQC
========================================================================================
*/

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMÈTRES DU PIPELINE
// ─────────────────────────────────────────────────────────────────────────────
params.fastq_dir    = null          // Dossier contenant les FASTQ (obligatoire)
params.outdir       = "output"      // Dossier de sortie des rapports HTML
params.resultsdir   = "results"     // Dossier de sortie des fichiers intermédiaires
params.min_length   = 200           // Longueur minimale des reads (NanoFilt)
params.min_quality  = 10            // Qualité minimale des reads Q-score (NanoFilt)
params.genome_size  = "5m"          // Taille estimée du génome pour Flye (ex: 5m = 5 Mb)
params.medaka_model = "r1041_e82_400bps_sup_v5.2.0"
                                    // Modèle Medaka : r1041 = R10.4.1 | e82 = Kit 14 | sup = SUP
params.bootstrap    = true          // Active le calcul des valeurs de bootstrap IQ-TREE
params.kraken_db    = "/data/kraken2_db"
                                    // Base de données Kraken2 (PlusPF-8 recommandée)
params.checkm2_db   = "/data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
                                    // Base de données CheckM2 (fichier .dmnd)

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// Chaque module = un outil = un fichier .nf dans le dossier modules/.
// ─────────────────────────────────────────────────────────────────────────────
include { NANOFILT }     from './modules/nanofilt.nf'     // Filtrage qualité des reads
include { NANOSTAT }     from './modules/nanostat.nf'     // Statistiques QC des reads
include { KRAKEN2 }      from './modules/kraken2.nf'      // Identification taxonomique
include { FLYE }         from './modules/flye.nf'         // Assemblage de novo
include { MEDAKA }       from './modules/medaka.nf'       // Polissage des assemblages
include { QUALIMAP }     from './modules/qualimap.nf'     // QC du mapping BAM
include { QUAST }        from './modules/quast.nf'        // QC structurel des assemblages
include { MLST }         from './modules/mlst.nf'         // Typage MLST
include { CHECKM2 }      from './modules/checkm2.nf'      // Complétude/contamination
include { KSNP4 }        from './modules/ksnp4.nf'        // SNP calling k-mer multi-souches
include { IQTREE }       from './modules/iqtree.nf'       // Arbre phylogénétique ML
include { PHYLO_REPORT } from './modules/phylo_report.nf' // Rapport HTML interactif
include { MULTIQC }      from './modules/multiqc.nf'      // Rapport QC agrégé

// ─────────────────────────────────────────────────────────────────────────────
// BANNIÈRE DE DÉMARRAGE
// Affiche les paramètres utilisés pour traçabilité dans les logs.
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
  Base Kraken2    : ${params.kraken_db}
  Base CheckM2    : ${params.checkm2_db}
──────────────────────────────────────────────────────────
""".stripIndent()

// ─────────────────────────────────────────────────────────────────────────────
// VALIDATION DES PARAMÈTRES OBLIGATOIRES
// On arrête le pipeline proprement si un paramètre requis est absent.
// ─────────────────────────────────────────────────────────────────────────────
if (!params.fastq_dir) {
    error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
}
if (!params.kraken_db) {
    error "ERREUR : --kraken_db est obligatoire.\nUsage : nextflow run main.nf --kraken_db /chemin/vers/db"
}
if (!params.checkm2_db) {
    error "ERREUR : --checkm2_db est obligatoire.\nUsage : nextflow run main.nf --checkm2_db /chemin/vers/uniref100.KO.1.dmnd"
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Création du channel d'entrée ──────────────────────────────────────────
    // LEÇON : fromPath() crée un channel à partir de fichiers sur le disque.
    // map() transforme chaque fichier en tuple [sample_id, fichier].
    // Le sample_id est extrait du nom de fichier en supprimant les extensions.
    // checkIfExists:true plante proprement si aucun fichier ne correspond.
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

    // ── Filtrage et QC des reads ───────────────────────────────────────────────
    // NanoFilt supprime les reads trop courts ou de mauvaise qualité.
    // NanoStat génère un rapport statistique par échantillon pour MultiQC.
    NANOFILT(ch_fastq)
    NANOSTAT(NANOFILT.out.reads)

    // ── Identification taxonomique ─────────────────────────────────────────────
    // Kraken2 classifie les reads contre la base PlusPF-8 pour confirmer
    // l'identité de l'espèce et détecter les contaminations (>1% des reads).
    // LEÇON : .first() transforme le channel en value channel réutilisable —
    // la référence est partagée entre toutes les souches sans être consommée.
    ch_kraken_db = Channel.fromPath(params.kraken_db, checkIfExists: true).first()
    KRAKEN2(NANOFILT.out.reads, ch_kraken_db)

    // ── Assemblage de novo ─────────────────────────────────────────────────────
    // Flye est optimisé pour les reads longs avec taux d'erreur élevé.
    FLYE(NANOFILT.out.reads)

    // ── Polissage des assemblages ──────────────────────────────────────────────
    // Medaka corrige les erreurs résiduelles via un modèle de réseau de neurones.
    // LEÇON : .join() garantit que chaque souche reçoit SES propres reads
    // et SON propre assemblage, en les associant par sample_id.
    ch_medaka_input = NANOFILT.out.reads.join(FLYE.out.assembly)
    MEDAKA(ch_medaka_input)

    // ── Contrôle qualité des assemblages ──────────────────────────────────────
    // Trois outils complémentaires lancés en parallèle sur les assemblages :
    //   - Qualimap : qualité du mapping BAM (couverture, profondeur)
    //   - QUAST    : qualité structurelle (N50, nb contigs, taille)
    //   - MLST     : typage séquence type (schéma PubMLST auto-détecté)
    //   - CheckM2  : complétude et contamination biologiques
    QUALIMAP(MEDAKA.out.bam)
    QUAST(MEDAKA.out.assembly)
    MLST(MEDAKA.out.assembly)
    ch_checkm2_db = Channel.fromPath(params.checkm2_db, checkIfExists: true).first()
    CHECKM2(MEDAKA.out.assembly, ch_checkm2_db)

    // ── SNP calling multi-souches ──────────────────────────────────────────────
    // kSNP4 compare tous les assemblages simultanément via une approche k-mer
    // (sans alignement global) — robuste aux réarrangements génomiques.
    // Produit deux alignements : tous les SNPs et SNPs core (présents dans
    // toutes les souches).
    // LEÇON : .map() extrait les FASTA sans le sample_id (kSNP4 le déduit
    // du nom de fichier). .collect() attend que toutes les souches soient
    // assemblées avant de lancer kSNP4.
    ch_all_assemblies = MEDAKA.out.assembly
        .map { sample_id, fasta -> fasta }
        .collect()

    KSNP4(ch_all_assemblies)

    // ── Arbre phylogénétique ───────────────────────────────────────────────────
    // IQ-TREE construit l'arbre par maximum de vraisemblance (GTR+G+ASC).
    // ASC = ascertainment bias correction, obligatoire pour un alignement
    // de SNPs uniquement (sites invariants exclus par kSNP4).
    IQTREE(KSNP4.out.snp_alignment)

    // ── Collecte des fichiers QC pour PHYLO_REPORT ────────────────────────────
    // On collecte tous les fichiers de chaque outil QC en une liste.
    // Nextflow les stage tous dans le même workdir de PHYLO_REPORT,
    // et le script Python lit ./ pour trouver les fichiers de chaque souche.
    ch_kraken_files  = KRAKEN2.out.report
        .map { sample_id, report -> report }
        .collect()

    ch_mlst_files    = MLST.out.mlst
        .map { sample_id, tsv -> tsv }
        .collect()

    ch_qualimap_dirs = QUALIMAP.out.results
        .collect()

    ch_checkm2_files = CHECKM2.out.report
        .map { sample_id, tsv -> tsv }
        .collect()

    // ── Rapports phylogénétiques interactifs ───────────────────────────────────
    // Deux rapports sont générés en parallèle :
    //   - report_all_snps  : basé sur tous les SNPs (vision globale)
    //   - report_core_snps : basé sur les SNPs core uniquement (plus conservateur,
    //                        positions présentes dans toutes les souches)
    // LEÇON : .combine() associe chaque FASTA avec l'arbre IQ-TREE.
    // .mix() fusionne les deux channels pour lancer PHYLO_REPORT deux fois
    // en parallèle.
    ch_all_snps = KSNP4.out.snp_alignment
        .map { fasta -> tuple("report_all_snps", fasta) }
        .combine(IQTREE.out.tree)

    ch_core_snps = KSNP4.out.core_snp_matrix
        .map { fasta -> tuple("report_core_snps", fasta) }
        .combine(IQTREE.out.tree)

    PHYLO_REPORT(
        ch_all_snps.mix(ch_core_snps),
        ch_kraken_files,
        ch_mlst_files,
        ch_qualimap_dirs,
        ch_checkm2_files
    )

    // ── Rapport MultiQC ────────────────────────────────────────────────────────
    // MultiQC agrège les rapports de tous les outils QC en un seul fichier HTML.
    // LEÇON : .mix() fusionne plusieurs channels en un seul flux.
    // .collect() attend que tous les fichiers soient disponibles avant de
    // lancer MultiQC.
    ch_multiqc_files = NANOSTAT.out.stats
        .mix(FLYE.out.stats)
        .mix(QUALIMAP.out.results)
        .mix(QUAST.out.results)
        .mix(CHECKM2.out.report.map { sample_id, report -> report })
        .mix(MLST.out.mlst.map { sample_id, tsv -> tsv })
        .collect()

    MULTIQC(ch_multiqc_files)
}

// ─────────────────────────────────────────────────────────────────────────────
// RÉSUMÉ DE FIN DE PIPELINE
// S'exécute automatiquement à la fin, succès ou échec.
// ─────────────────────────────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
╔══════════════════════════════════════════════════════════╗
║                  PIPELINE TERMINÉ !                      ║
╚══════════════════════════════════════════════════════════╝
  Statut    : ${workflow.success ? '✅ Succès' : '❌ Échec'}
  Durée     : ${workflow.duration}
  Rapports  : ${params.outdir}/
  Résultats : ${params.resultsdir}/
──────────────────────────────────────────────────────────
""".stripIndent()
}
