#!/usr/bin/env nextflow
/*
========================================================================================
    PIPELINE NANOPORE - ANALYSE DE SOUCHES BACTÉRIENNES v2.3
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    DEUX MODES D'EXÉCUTION :

    1. Mode de novo (sans --reference) :
       NANOFILT → NANOSTAT → FLYE → MEDAKA → QUALIMAP
                                           → KSNP4 → IQTREE → PHYLO_REPORT
                                           → MULTIQC

    2. Mode référence (avec --reference) :
       NANOFILT → NANOSTAT → MINIMAP2 → QUALIMAP
                                      → CLAIR3 → BCFTOOLS_MERGE → VCF_TO_FASTA
                                                                 → IQTREE → PHYLO_REPORT
                                      → MULTIQC
========================================================================================
*/

// LEÇON : DSL2 est la version moderne de Nextflow.
// Elle permet de réutiliser des modules et d'avoir un code plus propre.
// Activé par défaut depuis Nextflow 22.03, mais on le garde pour la lisibilité.
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMÈTRES DU PIPELINE
// L'utilisateur peut les modifier en ligne de commande :
//   nextflow run main.nf --fastq_dir /mon/dossier --reference genome.fasta
// EPI2ME les expose dans son interface via nextflow_schema.json.
// ─────────────────────────────────────────────────────────────────────────────
params.fastq_dir    = null          // Dossier contenant les FASTQ (obligatoire)
params.reference    = null          // Génome de référence FASTA (optionnel)
                                    // Si fourni → mode référence (Minimap2 + Clair3)
                                    // Si absent → mode de novo (Flye + Medaka + kSNP4)
params.outdir       = "output"      // Dossier de sortie
params.min_length   = 200           // Longueur minimale des reads (NanoFilt)
params.min_quality  = 10            // Qualité minimale des reads Q-score (NanoFilt)
params.genome_size  = "5m"          // Taille estimée du génome pour Flye (ex: 5m = 5 Mb)
                                    // Utilisé uniquement en mode de novo
params.medaka_model = "r1041_e82_400bps_sup_v5.2.0"
                                    // Modèle Medaka : doit correspondre à ta flowcell
                                    // et au mode de basecalling utilisé dans MinKNOW
                                    // r1041 = R10.4.1 | e82 = Kit 14 | sup = SUP
                                    // Utilisé uniquement en mode de novo
params.clair3_model = "r1041_e82_400bps_sup_v500"
                                    // Modèle Clair3 : doit correspondre à ta chimie ONT
                                    // Convention : v500 = Dorado SUP 5.0.0
                                    // Modèles disponibles dans /opt/models/ du container
                                    // Utilisé uniquement en mode référence
params.bootstrap    = true          // Active le calcul des valeurs de bootstrap IQ-TREE
                                    // Désactiver (false) pour accélérer les tests

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// Chaque module = un outil = un fichier .nf dans le dossier modules/.
// Cette séparation rend le code modulaire et facilite la maintenance.
// ─────────────────────────────────────────────────────────────────────────────

// Modules communs aux deux modes
include { NANOFILT }       from './modules/nanofilt.nf'       // Filtrage qualité des reads
include { NANOSTAT }       from './modules/nanostat.nf'       // Statistiques QC des reads
include { QUALIMAP }       from './modules/qualimap.nf'       // QC du mapping BAM
include { IQTREE }         from './modules/iqtree.nf'         // Construction de l'arbre phylogénétique
include { MULTIQC }        from './modules/multiqc.nf'        // Rapport QC agrégé
include { PHYLO_REPORT }   from './modules/phylo_report.nf'   // Rapport HTML interactif

// Modules mode de novo uniquement
include { FLYE }           from './modules/flye.nf'           // Assemblage de novo
include { MEDAKA }         from './modules/medaka.nf'         // Polissage des assemblages
include { KSNP4 }          from './modules/ksnp4.nf'          // SNP calling sans référence (k-mer)

// Modules mode référence uniquement
include { MINIMAP2 }       from './modules/minimap2.nf'       // Mapping long reads → référence
include { CLAIR3 }         from './modules/clair3.nf'         // Variant calling deep learning
include { BCFTOOLS_MERGE } from './modules/bcftools_merge.nf' // Fusion des VCF multi-souches
include { VCF_TO_FASTA }   from './modules/vcf_to_fasta.nf'  // Conversion VCF → alignement FASTA

// ─────────────────────────────────────────────────────────────────────────────
// BANNIÈRE DE DÉMARRAGE
// Affiche les paramètres utilisés pour traçabilité dans les logs.
// ─────────────────────────────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════╗
║     PIPELINE SNP & PHYLOGÉNIE - NANOPORE MINION  v2.3    ║
╚══════════════════════════════════════════════════════════╝
  Dossier FASTQ   : ${params.fastq_dir}
  Référence       : ${params.reference ?: 'Aucune — mode assemblage de novo'}
  Dossier sortie  : ${params.outdir}
  Qualité min.    : ${params.min_quality}
  Longueur min.   : ${params.min_length} bp
  Taille génome   : ${params.genome_size}
  Modèle Medaka   : ${params.medaka_model}
  Modèle Clair3   : ${params.clair3_model}
──────────────────────────────────────────────────────────
""".stripIndent()

// ─────────────────────────────────────────────────────────────────────────────
// VALIDATION DES PARAMÈTRES OBLIGATOIRES
// On arrête le pipeline proprement si l'utilisateur oublie --fastq_dir.
// ─────────────────────────────────────────────────────────────────────────────
if (!params.fastq_dir) {
    error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
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

    // ── Étapes communes aux deux modes ────────────────────────────────────────
    // Ces étapes s'exécutent toujours, quelle que soit la présence d'une référence.

    // Filtrage qualité : supprime les reads trop courts ou de mauvaise qualité
    NANOFILT(ch_fastq)

    // Statistiques : génère un rapport qualité par échantillon
    NANOSTAT(NANOFILT.out.reads)

    // ════════════════════════════════════════════════════════════════════════════
    // MODE RÉFÉRENCE (params.reference fourni)
    // Approche : mapping des reads sur la référence → variant calling par souche
    // → fusion des VCF → alignement FASTA multi-souches → phylogénie
    //
    // Avantages vs mode de novo :
    //   - Plus rapide (pas d'assemblage)
    //   - Meilleure précision SNP (Clair3 deep learning vs kSNP4 k-mer)
    //   - Variants exprimés par rapport à une référence connue
    // ════════════════════════════════════════════════════════════════════════════
    if (params.reference) {

        // Chargement de la référence comme channel singleton
        // LEÇON : un channel singleton émet sa valeur à chaque fois qu'un
        // process en a besoin — idéal pour une référence partagée entre toutes
        // les souches sans avoir à la dupliquer dans le code.
        ch_reference = Channel.fromPath(params.reference, checkIfExists: true)

        // Mapping : aligne les reads filtrés sur la référence avec Minimap2
        // preset map-ont = optimisé pour les reads longs Oxford Nanopore
        // Produit un BAM trié et indexé par souche
        MINIMAP2(NANOFILT.out.reads, ch_reference)

        // QC du mapping : vérifie la couverture et la qualité de l'alignement
        // Les rapports HTML sont nommés par sample_id pour identification dans EPI2ME
        QUALIMAP(MINIMAP2.out.bam)

        // Variant calling : détecte les SNPs par souche via réseau de neurones
        // Clair3 utilise une approche double : pileup (rapide) + full-alignment (précis)
        // Le modèle doit correspondre à la chimie ONT et au mode de basecalling
        CLAIR3(MINIMAP2.out.bam, ch_reference)

        // Collecte tous les VCF individuels pour les merger
        // LEÇON : .map() extrait uniquement le VCF du tuple (on n'a plus besoin
        // du sample_id à ce stade). .collect() attend que TOUTES les souches
        // soient terminées avant de passer à l'étape suivante.
        ch_all_vcf = CLAIR3.out.vcf
            .map { sample_id, vcf -> vcf }
            .collect()

        // Fusion : combine tous les VCF individuels en un VCF multi-sample
        // bcftools merge aligne les positions variables entre les souches
        // et gère les génotypes manquants (./.) pour les positions non appelées
        BCFTOOLS_MERGE(ch_all_vcf)

        // Conversion : transforme le VCF multi-sample en alignement FASTA
        // Chaque séquence = une souche, chaque position = un SNP variable
        // Ce FASTA est le format attendu par IQ-TREE et PHYLO_REPORT
        VCF_TO_FASTA(BCFTOOLS_MERGE.out.merged_vcf)

        // Arbre phylogénétique à partir de l'alignement SNP
        IQTREE(VCF_TO_FASTA.out.fasta)

        // Rapport interactif HTML : arbre + matrice de distances SNP
        // LEÇON : .combine() associe le FASTA avec l'arbre IQ-TREE
        // pour les passer ensemble à PHYLO_REPORT dans un seul tuple.
        ch_phylo = VCF_TO_FASTA.out.fasta
            .map { fasta -> tuple("report_reference_snps", fasta) }
            .combine(IQTREE.out.tree)

        PHYLO_REPORT(ch_phylo)

        // Rapport MultiQC : agrège les stats NanoStat + Qualimap
        // LEÇON : .mix() fusionne plusieurs channels en un seul flux.
        // .collect() attend tous les fichiers avant de lancer MultiQC.
        ch_multiqc_files = NANOSTAT.out.stats
            .mix(QUALIMAP.out.results)
            .collect()

        MULTIQC(ch_multiqc_files)

    // ════════════════════════════════════════════════════════════════════════════
    // MODE DE NOVO (pas de référence fournie)
    // Approche : assemblage de chaque souche → polissage → comparaison
    // multi-souches par k-mers (kSNP4) → phylogénie
    //
    // Avantages vs mode référence :
    //   - Pas besoin de référence (utile pour espèces peu caractérisées)
    //   - Capture les variants structuraux absents de la référence
    //   - Produit des assemblages complets réutilisables
    // ════════════════════════════════════════════════════════════════════════════
    } else {

        // Assemblage de novo : reconstruit le génome depuis les reads
        // Flye est optimisé pour les reads longs avec taux d'erreur élevé
        FLYE(NANOFILT.out.reads)

        // Polissage : corrige les erreurs résiduelles de l'assemblage
        // LEÇON : .join() combine les reads et l'assemblage pour chaque
        // souche sur la base du sample_id — garantit que chaque souche
        // reçoit bien SES propres reads et SON propre assemblage.
        ch_medaka_input = NANOFILT.out.reads.join(FLYE.out.assembly)
        MEDAKA(ch_medaka_input)

        // QC du mapping : Medaka produit un BAM d'alignement interne
        // qu'on utilise pour vérifier la couverture de l'assemblage
        QUALIMAP(MEDAKA.out.bam)

        // SNP calling multi-souches : compare tous les assemblages simultanément
        // kSNP4 utilise une approche k-mer (sans alignement global) :
        //   - Plus rapide que l'alignement pour de nombreuses souches
        //   - Robuste aux réarrangements génomiques
        //   - Produit un alignement SNP (tous SNPs) et un alignement SNP core
        // LEÇON : .map() extrait les FASTA sans le sample_id (kSNP4 le
        // déduit du nom de fichier). .collect() attend toutes les souches.
        ch_all_assemblies = MEDAKA.out.assembly
            .map { sample_id, fasta -> fasta }
            .collect()

        KSNP4(ch_all_assemblies)

        // Arbre phylogénétique à partir de l'alignement SNP kSNP4
        IQTREE(KSNP4.out.snp_alignment)

        // Deux rapports interactifs :
        //   - Tous les SNPs (snp_alignment.fasta) : vision globale
        //   - SNPs core uniquement (core_snp_matrix.fasta) : positions
        //     présentes dans TOUTES les souches, plus conservateur
        // LEÇON : .mix() fusionne les deux channels pour lancer
        // PHYLO_REPORT deux fois en parallèle (une fois par rapport).
        ch_all_snps = KSNP4.out.snp_alignment
            .map { fasta -> tuple("report_all_snps", fasta) }
            .combine(IQTREE.out.tree)

        ch_core_snps = KSNP4.out.core_snp_matrix
            .map { fasta -> tuple("report_core_snps", fasta) }
            .combine(IQTREE.out.tree)

        PHYLO_REPORT(ch_all_snps.mix(ch_core_snps))

        // Rapport MultiQC : agrège NanoStat + stats Flye + Qualimap
        ch_multiqc_files = NANOSTAT.out.stats
            .mix(FLYE.out.stats)
            .mix(QUALIMAP.out.results)
            .collect()

        MULTIQC(ch_multiqc_files)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RÉSUMÉ DE FIN DE PIPELINE
// S'exécute automatiquement à la fin, succès ou échec.
// ─────────────────────────────────────────────────────────────────────────────
workflow.onComplete {
    def mode = params.reference ? 'Référence' : 'De novo'
    log.info """
╔══════════════════════════════════════════════════════════╗
║                  PIPELINE TERMINÉ !                      ║
╚══════════════════════════════════════════════════════════╝
  Statut    : ${workflow.success ? '✅ Succès' : '❌ Échec'}
  Mode      : ${mode}
  Durée     : ${workflow.duration}
  Résultats : ${params.outdir}/
──────────────────────────────────────────────────────────
""".stripIndent()
}
