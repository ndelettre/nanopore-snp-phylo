#!/usr/bin/env nextflow
/*
========================================================================================
    PIPELINE NANOPORE - ANALYSE DE SOUCHES BACTÉRIENNES
    (version : manifest.version dans nextflow.config)
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    WORKFLOW :
       CHOPPER → NANOSTAT → KRAKEN2 → BRACKEN (identification espèce)
                          → FLYE → MEDAKA → QUALIMAP
                                          → QUAST
                                          → MLST
                                          → CHECKM2
                                          → KSNP4 → IQTREE ×2 → PHYLO_REPORT ×2
                                            (tous les SNPs / SNPs core)
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
params.min_length   = 1000          // Longueur minimale des reads (Chopper)
params.min_quality  = 10            // Qualité minimale des reads Q-score (Chopper)
params.genome_size  = "5m"          // Taille estimée du génome pour Flye (ex: 5m = 5 Mb)
params.medaka_model = "r1041_e82_400bps_bacterial_methylation"
                                    // Modèle Medaka : r1041 = R10.4.1 | e82 = Kit 14 | sup = SUP
params.bootstrap    = false         // Active le calcul des valeurs de bootstrap IQ-TREE (--bootstrap)
params.kraken_db    = "/data/kraken2_db"
                                    // Base de données Kraken2 (PlusPF-16 recommandée)
params.checkm2_db   = "/data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
                                    // Base de données CheckM2 (fichier .dmnd)
params.reference_fasta = null       // Optionnel : assemblage de référence (FASTA non compressé)

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// Chaque module = un outil = un fichier .nf dans le dossier modules/.
// ─────────────────────────────────────────────────────────────────────────────
include { CHOPPER }      from './modules/chopper.nf'      // Filtrage qualité des reads
include { NANOSTAT }     from './modules/nanostat.nf'     // Statistiques QC des reads
include { KRAKEN2 }      from './modules/kraken2.nf'      // Identification taxonomique
include { KRAKEN2_REFERENCE } from './modules/kraken2.nf' // Espèce de la référence
include { BRACKEN }      from './modules/bracken.nf'      // Réestimation au rang espèce
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
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // LEÇON : en syntaxe stricte (Nextflow ≥ 25), aucune instruction ne doit
    // se trouver hors d'un bloc workflow ou process. La bannière, la
    // validation des paramètres et le résumé de fin sont donc ici.

    // ── Bannière de démarrage ─────────────────────────────────────────────────
    // Affiche les paramètres utilisés pour traçabilité dans les logs.
    log.info """
╔══════════════════════════════════════════════════════════╗
║     PIPELINE SNP & PHYLOGÉNIE - NANOPORE MINION  v${workflow.manifest.version.padRight(7)}║
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
  Référence       : ${params.reference_fasta ?: 'aucune'}
──────────────────────────────────────────────────────────
""".stripIndent()

    // ── Validation des paramètres obligatoires ────────────────────────────────
    // On arrête le pipeline proprement si un paramètre requis est absent.
    if (!params.fastq_dir) {
        error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
    }
    if (!params.kraken_db) {
        error "ERREUR : --kraken_db est obligatoire.\nUsage : nextflow run main.nf --kraken_db /chemin/vers/db"
    }
    if (!params.checkm2_db) {
        error "ERREUR : --checkm2_db est obligatoire.\nUsage : nextflow run main.nf --checkm2_db /chemin/vers/uniref100.KO.1.dmnd"
    }

    // ── Résumé de fin de pipeline ─────────────────────────────────────────────
    // S'exécute automatiquement à la fin, succès ou échec.
    // LEÇON : dans le bloc workflow, `workflow` et `params` ne sont plus
    // visibles depuis la closure du handler → on les garde dans des
    // variables locales, que la closure capture.
    def run_info   = workflow
    def outdir     = params.outdir
    def resultsdir = params.resultsdir
    workflow.onComplete {
        log.info """
╔══════════════════════════════════════════════════════════╗
║                  PIPELINE TERMINÉ !                      ║
╚══════════════════════════════════════════════════════════╝
  Statut    : ${run_info.success ? '✅ Succès' : '❌ Échec'}
  Durée     : ${run_info.duration}
  Rapports  : ${outdir}/
  Résultats : ${resultsdir}/
──────────────────────────────────────────────────────────
""".stripIndent()
    }

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

    // ── Souche de référence optionnelle ───────────────────────────────────────
    // Assemblage déjà fait, fourni via --reference_fasta : il rejoint MLST et
    // kSNP4, mais pas le QC (Qualimap, QUAST, CheckM2). Kraken2 ne sert qu'à
    // afficher son espèce dans le rapport.
    // Nom de la souche = nom du fichier jusqu'au premier point.
    // LEÇON : sans référence, Channel.empty() → les .mix() plus bas n'ajoutent
    // rien et le pipeline se comporte exactement comme avant.
    if (params.reference_fasta) {
        if (params.reference_fasta.toString() =~ /\.gz$/) {
            error "ERREUR : --reference_fasta doit être un FASTA non compressé (kSNP4 ne lit pas le .gz)."
        }
        ch_reference = Channel
            .fromPath(params.reference_fasta, checkIfExists: true)
            .map { fasta -> tuple(fasta.getSimpleName(), fasta) }
    } else {
        ch_reference = Channel.empty()
    }

    // ── Filtrage et QC des reads ───────────────────────────────────────────────
    // Chopper supprime les reads trop courts ou de mauvaise qualité.
    // NanoStat génère un rapport statistique par échantillon pour MultiQC.
    CHOPPER(ch_fastq)
    NANOSTAT(CHOPPER.out.reads)

    // ── Identification taxonomique ─────────────────────────────────────────────
    // Kraken2 classifie les reads contre la base PlusPF-16 pour confirmer
    // l'identité de l'espèce et détecter les contaminations (>1% des reads).
    // LEÇON : .first() transforme le channel en value channel réutilisable —
    // la référence est partagée entre toutes les souches sans être consommée.
    ch_kraken_db = Channel.fromPath(params.kraken_db, checkIfExists: true).first()
    KRAKEN2(CHOPPER.out.reads, ch_kraken_db)

    // Bracken redistribue au rang espèce les reads que Kraken2 a laissés au
    // rang genre (reads communs à plusieurs espèces proches).
    BRACKEN(KRAKEN2.out.report, ch_kraken_db)

    // Espèce de la référence : classification de ses contigs.
    KRAKEN2_REFERENCE(ch_reference, ch_kraken_db)

    // ── Assemblage de novo ─────────────────────────────────────────────────────
    // Flye est optimisé pour les reads longs avec taux d'erreur élevé.
    FLYE(CHOPPER.out.reads)

    // ── Polissage des assemblages ──────────────────────────────────────────────
    // Medaka corrige les erreurs résiduelles via un modèle de réseau de neurones.
    // LEÇON : .join() garantit que chaque souche reçoit SES propres reads
    // et SON propre assemblage, en les associant par sample_id.
    ch_medaka_input = CHOPPER.out.reads.join(FLYE.out.assembly)
    MEDAKA(ch_medaka_input)

    // ── Contrôle qualité des assemblages ──────────────────────────────────────
    // Quatre outils complémentaires lancés en parallèle sur les assemblages :
    //   - Qualimap : qualité du mapping BAM (couverture, profondeur)
    //   - QUAST    : qualité structurelle (N50, nb contigs, taille)
    //   - MLST     : typage séquence type (schéma PubMLST auto-détecté)
    //   - CheckM2  : complétude et contamination biologiques
    QUALIMAP(MEDAKA.out.bam)
    QUAST(MEDAKA.out.assembly)
    MLST(MEDAKA.out.assembly.mix(ch_reference))
    ch_checkm2_db = Channel.fromPath(params.checkm2_db, checkIfExists: true).first()
    CHECKM2(MEDAKA.out.assembly, ch_checkm2_db)

    // ── SNP calling multi-souches ──────────────────────────────────────────────
    // kSNP4 compare tous les assemblages simultanément via une approche k-mer
    // (sans alignement global) — robuste aux réarrangements génomiques.
    // Produit deux alignements : tous les SNPs (SNPs_all_matrix.fasta) et
    // SNPs core (présents dans toutes les souches).
    // LEÇON : .map() extrait les FASTA sans le sample_id (kSNP4 le déduit
    // du nom de fichier). .collect() attend que toutes les souches soient
    // assemblées avant de lancer kSNP4.
    ch_all_assemblies = MEDAKA.out.assembly
        .mix(ch_reference)
        .map { sample_id, fasta -> fasta }
        .collect()

    KSNP4(ch_all_assemblies)

    // ── Arbres phylogénétiques ─────────────────────────────────────────────────
    // IQ-TREE construit un arbre par maximum de vraisemblance (GTR+G+ASC)
    // pour chaque alignement. ASC = ascertainment bias correction, obligatoire
    // pour un alignement de SNPs uniquement (sites invariants exclus par kSNP4).
    // Chaque alignement porte le nom du rapport qu'il alimente :
    //   - report_all_snps  : tous les SNPs (distances fiables même avec une
    //                        souche éloignée, positions manquantes ignorées)
    //   - report_core_snps : SNPs core uniquement (plus conservateur,
    //                        positions présentes dans toutes les souches)
    ch_alignments = KSNP4.out.all_snp_alignment
        .map { fasta -> tuple("report_all_snps", fasta) }
        .mix(
            KSNP4.out.core_snp_alignment
                .map { fasta -> tuple("report_core_snps", fasta) }
        )

    IQTREE(ch_alignments)

    // ── Collecte des fichiers QC pour PHYLO_REPORT ────────────────────────────
    // On collecte tous les fichiers de chaque outil QC en une liste.
    // Nextflow les stage tous dans le même workdir de PHYLO_REPORT,
    // et le script Python lit ./ pour trouver les fichiers de chaque souche.
    ch_kraken_files  = BRACKEN.out.abundance
        .map { sample_id, tsv -> tsv }
        .collect()

    // LEÇON : sans référence le channel est vide, et .collect() n'émettrait
    // rien → PHYLO_REPORT ne serait jamais lancé. .ifEmpty([]) fournit une
    // liste vide à la place.
    ch_ref_kraken_files = KRAKEN2_REFERENCE.out.classification
        .map { sample_id, report, out -> [report, out] }
        .collect()
        .ifEmpty([])

    ch_mlst_files    = MLST.out.mlst
        .map { sample_id, tsv -> tsv }
        .collect()

    ch_qualimap_dirs = QUALIMAP.out.results
        .collect()

    ch_checkm2_files = CHECKM2.out.report
        .map { sample_id, tsv -> tsv }
        .collect()

    // ── Rapports phylogénétiques interactifs ───────────────────────────────────
    // Un rapport par alignement, chacun avec son propre arbre.
    // LEÇON : .join() associe chaque alignement à SON arbre par le nom du
    // rapport → [report_name, fasta, treefile].
    ch_reports = ch_alignments.join(IQTREE.out.tree)

    PHYLO_REPORT(
        ch_reports,
        ch_kraken_files,
        ch_ref_kraken_files,
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
