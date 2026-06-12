#!/usr/bin/env nextflow

/*
========================================================================================
    PIPELINE NANOPORE - ANALYSE DE SOUCHES BACTÉRIENNES v2.2 (corrigé)
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    LEÇON : Ce fichier est le point d'entrée du pipeline.
    Il importe les modules et définit le workflow principal.
========================================================================================
*/

// LEÇON : DSL2 est la version moderne de Nextflow.
// Elle permet de réutiliser des modules et d'avoir un code plus propre.
// NOTE : depuis Nextflow 22.03, DSL2 est activé par défaut, cette ligne est
// facultative mais on la garde pour la compatibilité avec d'anciennes versions.
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// LEÇON : Les "params" sont les paramètres du pipeline.
// L'utilisateur peut les modifier en ligne de commande :
//   nextflow run main.nf --fastq_dir /mon/dossier --reference genome.fasta
// ─────────────────────────────────────────────────────────────────────────────
params.fastq_dir    = null          // Dossier contenant les FASTQ (obligatoire)
params.reference    = null          // Génome de référence FASTA (optionnel)
params.outdir       = "results"     // Dossier de sortie
params.min_length   = 1000          // Longueur minimale des reads (NanoFilt)
params.min_quality  = 10            // Qualité minimale des reads (NanoFilt)
params.genome_size  = "5m"          // Taille estimée du génome pour Flye (ex: 5m = 5 Mb)
params.medaka_model = "r1041_e82_400bps_sup_v5.0.0"  // Modèle Medaka selon ta flowcell
params.bootstrap = true

// ─────────────────────────────────────────────────────────────────────────────
// LEÇON : On importe les modules depuis le dossier modules/.
// Chaque module = un outil = un fichier .nf séparé.
// Cela rend le code modulaire et réutilisable.
//
// CORRECTION v2.1 : on ajoute l'extension .nf explicite à chaque include.
// Nextflow accepte les deux formes, mais la forme avec .nf est plus claire
// et évite toute ambiguïté si un sous-dossier portait le même nom qu'un module.
// ─────────────────────────────────────────────────────────────────────────────
include { NANOFILT }  from './modules/nanofilt.nf'
include { NANOSTAT }  from './modules/nanostat.nf'
include { FLYE }      from './modules/flye.nf'
include { MEDAKA }    from './modules/medaka.nf'
include { KSNP4 }     from './modules/ksnp4.nf'
include { SNIPPY }    from './modules/snippy.nf'
include { IQTREE }    from './modules/iqtree.nf'
include { MULTIQC }   from './modules/multiqc.nf'
include { QUALIMAP }  from './modules/qualimap.nf'

// ─────────────────────────────────────────────────────────────────────────────
// LEÇON : La fonction "log.info" affiche un message au démarrage du pipeline.
// C'est une bonne pratique pour documenter les paramètres utilisés.
// ─────────────────────────────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════╗
║     PIPELINE SNP & PHYLOGÉNIE - NANOPORE MINION  v2.2    ║
╚══════════════════════════════════════════════════════════╝
  Dossier FASTQ   : ${params.fastq_dir}
  Référence       : ${params.reference ?: 'Aucune (mode de novo)'}
  Dossier sortie  : ${params.outdir}
  Qualité min.    : ${params.min_quality}
  Longueur min.   : ${params.min_length} bp
  Taille génome   : ${params.genome_size}
  Modèle Medaka   : ${params.medaka_model}
──────────────────────────────────────────────────────────
""".stripIndent()

// ─────────────────────────────────────────────────────────────────────────────
// LEÇON : Validation des paramètres obligatoires.
// On arrête le pipeline proprement si l'utilisateur oublie --fastq_dir.
// ─────────────────────────────────────────────────────────────────────────────
if (!params.fastq_dir) {
    error """
    ERREUR : Le paramètre --fastq_dir est obligatoire.
    Usage : nextflow run main.nf --fastq_dir /chemin/vers/fastq
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : Un "Channel" est un flux de données entre les process.
    // fromPath() crée un channel à partir de fichiers sur le disque.
    // map() transforme chaque élément : ici on crée une paire [nom_echantillon, fichier]
    // Cette paire s'appelle un "tuple" et c'est la convention standard dans Nextflow.
    //
    // Exemple : si tu as sample1.fastq.gz, le tuple sera :
    //   ["sample1", /chemin/vers/sample1.fastq.gz]
    //
    // CORRECTION v2.1 : on ajoute checkIfExists:true pour que Nextflow plante
    // proprement si le dossier n'existe pas, au lieu de continuer avec un
    // channel silencieusement vide. Cela remplace avantageusement l'usage
    // incorrect de .ifEmpty{} qu'on avait auparavant.
    // ─────────────────────────────────────────────────────────────────────────
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

    // ─────────────────────────────────────────────────────────────────────────
    // CORRECTION v2.1 (bug majeur) : la version précédente faisait
    //     ch_fastq.ifEmpty { error "..." }
    // C'était incorrect pour deux raisons :
    //   1) .ifEmpty() est un OPÉRATEUR qui renvoie un nouveau channel avec
    //      la valeur par défaut s'il est vide, il ne s'exécute pas comme
    //      une vérification "sur place". La closure {error ...} n'aurait
    //      jamais été évaluée au bon moment.
    //   2) En DSL2, un channel queue ne peut être consommé qu'UNE SEULE FOIS.
    //      Faire ch_fastq.ifEmpty{...} puis NANOFILT(ch_fastq) est interdit.
    //
    // Solution : on utilise checkIfExists:true dans fromPath (voir ci-dessus),
    // ce qui déclenche une erreur Nextflow si le pattern ne matche rien.
    // Plus simple et conforme aux idiomes DSL2.
    // ─────────────────────────────────────────────────────────────────────────

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 1 : Filtrage qualité avec NanoFilt
    // Input  : [sample_id, fastq_brut]
    // Output : [sample_id, fastq_filtré]
    // ─────────────────────────────────────────────────────────────────────────
    NANOFILT(ch_fastq)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 2 : Statistiques qualité avec NanoStat
    // Input  : [sample_id, fastq_filtré]
    // Output : rapport texte par échantillon
    //
    // CORRECTION v2.1 (commentaire factuel) : le commentaire précédent
    // affirmait que NanoStat est "nativement parsé par MultiQC". C'est FAUX :
    // MultiQC n'a pas de module natif pour NanoStat. Les fichiers seront
    // affichés via la détection "custom content" de MultiQC mais pas
    // intégrés en graphes natifs. Pour des graphes natifs, préférer
    // NanoPlot + un module MultiQC custom, ou pycoQC.
    // ─────────────────────────────────────────────────────────────────────────
    NANOSTAT(NANOFILT.out.reads)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 3 : Assemblage de novo avec Flye
    // Input  : [sample_id, fastq_filtré]
    // Output : [sample_id, assembly.fasta]
    // ─────────────────────────────────────────────────────────────────────────
    FLYE(NANOFILT.out.reads)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 4 : Polissage avec Medaka
    // Input  : [sample_id, fastq_filtré, assembly.fasta]
    // Output : [sample_id, assembly_polished.fasta]
    //
    // LEÇON : .join() combine deux channels sur la base de la clé (1er élément
    // du tuple par défaut, ici sample_id). C'est essentiel pour que chaque
    // échantillon reçoive SON propre assembly, et non celui d'un autre.
    // ─────────────────────────────────────────────────────────────────────────
    ch_medaka_input = NANOFILT.out.reads.join(FLYE.out.assembly)
    MEDAKA(ch_medaka_input)

    QUALIMAP(MEDAKA.out.bam)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 5 : kSNP4 sur tous les assemblages polishés
    // kSNP4 prend TOUS les génomes en même temps (pas échantillon par échantillon).
    //
    // LEÇON : .collect() attend que TOUS les échantillons soient finis,
    // puis les regroupe dans une seule liste. Indispensable pour kSNP4.
    //
    // LEÇON : .map { sample_id, fasta -> fasta } extrait uniquement le fichier
    // FASTA du tuple, car kSNP4 n'a pas besoin du sample_id (il le déduit
    // du nom de fichier, ou on le fournit dans le in_list).
    // ─────────────────────────────────────────────────────────────────────────
    ch_all_assemblies = MEDAKA.out.assembly
        .map { sample_id, fasta -> fasta }
        .collect()

    KSNP4(ch_all_assemblies)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 6 (optionnelle) : Si une référence est fournie → Snippy
    // LEÇON : Le bloc "if" conditionnel permet d'activer des branches du
    // pipeline selon les paramètres de l'utilisateur.
    //
    // LEÇON : .combine() crée toutes les combinaisons entre deux channels.
    // Ici : chaque assemblage polished × la référence unique = un tuple par souche.
    //
    // CORRECTION v2.1 : on ajoute checkIfExists:true sur la référence également.
    // ─────────────────────────────────────────────────────────────────────────
    if (params.reference) {
        ch_reference    = Channel.fromPath(params.reference, checkIfExists: true)
        ch_snippy_input = MEDAKA.out.assembly.combine(ch_reference)
        SNIPPY(ch_snippy_input)
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 7 : Arbre phylogénétique avec IQ-TREE
    // Utilise l'alignement SNP au format FASTA produit par kSNP4
    // ─────────────────────────────────────────────────────────────────────────
    IQTREE(KSNP4.out.snp_alignment)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 8 : Rapport MultiQC
    // LEÇON : .mix() fusionne plusieurs channels en un seul (union).
    // MultiQC lit tous les fichiers de stats et produit un rapport HTML agrégé.
    // ─────────────────────────────────────────────────────────────────────────
    ch_multiqc_files = NANOSTAT.out.stats
        .mix(FLYE.out.stats)
        .mix(QUALIMAP.out.results)
        .collect()

    MULTIQC(ch_multiqc_files)
}

// ─────────────────────────────────────────────────────────────────────────────
// LEÇON : Le bloc "workflow.onComplete" s'exécute à la fin du pipeline.
// C'est pratique pour afficher un résumé ou envoyer une notification.
// ─────────────────────────────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
╔══════════════════════════════════════════════════════════╗
║                  PIPELINE TERMINÉ !                      ║
╚══════════════════════════════════════════════════════════╝
  Statut    : ${workflow.success ? '✅ Succès' : '❌ Échec'}
  Durée     : ${workflow.duration}
  Résultats : ${params.outdir}/
──────────────────────────────────────────────────────────
  Fichiers clés :
    Matrice SNP    → ${params.outdir}/ksnp4/SNPs_all
    Arbre IQ-TREE  → ${params.outdir}/iqtree/phylo_tree.treefile
    Rapport QC     → ${params.outdir}/multiqc/multiqc_report.html
──────────────────────────────────────────────────────────
""".stripIndent()
}
