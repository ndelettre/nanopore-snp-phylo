/*
========================================================================================
    MODULE : kSNP4 — Détection de SNPs et phylogénie sans alignement global
========================================================================================
    kSNP4 identifie les SNPs à partir de k-mers (fragments de séquences de taille k).
    Il ne nécessite PAS de génome de référence et fonctionne sur des génomes complets.

    Sorties utilisées (dossier -outdir) :
    - SNPs_all_matrix.fasta   → alignement FASTA de TOUS les SNPs ; une souche
                                qui n'a pas le locus porte "-" (donnée manquante)
    - core_SNPs_matrix.fasta  → alignement FASTA des SNPs core (loci présents
                                dans TOUTES les souches)

    LEÇON : les distances SNP entre deux souches se calculent sur l'alignement
    de tous les SNPs, en ignorant les positions manquantes. Sur l'alignement
    core, un génome éloigné (ex. une référence) réduit le génome commun et
    fait disparaître des SNPs réels entre les autres souches.

    Les arbres sont construits par IQ-TREE : on ne demande pas à kSNP4 ses
    propres arbres ML/NJ (-ML, -NJ), qui ne seraient pas utilisés.
========================================================================================
*/

process KSNP4 {

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : Pas de tag ici car kSNP4 traite TOUS les échantillons ensemble.
    // Il n'y a donc pas de "sample_id" individuel.
    // ─────────────────────────────────────────────────────────────────────────
    label 'process_high'

    publishDir "${params.resultsdir}/ksnp4", mode: 'copy'

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : L'input est une liste de fichiers FASTA (path).
    // Cette liste a été créée dans main.nf par :
    //   MEDAKA.out.assembly
    //     .map { sample_id, fasta -> fasta }   ← extrait uniquement le FASTA
    //     .collect()                            ← attend tous les échantillons
    //
    // On n'utilise pas de tuple val/path ici car .collect() sur des tuples
    // produirait une liste plate que Nextflow ne saurait pas interpréter.
    // ─────────────────────────────────────────────────────────────────────────
    input:
    path fastas

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : optional:true sur un output signifie que Nextflow ne plantera
    // pas si ce fichier n'est pas produit. Sans aucun SNP core (souches trop
    // éloignées), le rapport core n'est simplement pas généré.
    // ─────────────────────────────────────────────────────────────────────────
    output:
    path "SNPs_all_matrix.fasta",                     emit: all_snp_alignment
    path "core_SNPs_matrix.fasta", optional: true,    emit: core_snp_alignment
    path "ksnp4_results/",                            emit: all_results

    script:
    """
    set -euo pipefail

    # ─────────────────────────────────────────────────────────────────────────
    # kSNP4 requiert un fichier "in_list" au format :
    #   /chemin/absolu/vers/genome.fasta<TAB>NomEchantillon
    # On génère ce fichier dynamiquement à partir des fichiers reçus.
    # Nom = nom du fichier jusqu'au premier point (même règle que getSimpleName()
    # dans main.nf), puis sed supprime le suffixe "_polished" des assemblages
    # Medaka. Une référence ajoutée (.fasta, .fna, .fa…) obtient ainsi le même
    # nom que son rapport MLST, ce qui permet la jointure dans PHYLO_REPORT.
    # ─────────────────────────────────────────────────────────────────────────
    > genome_list.txt
    for fasta in ${fastas}; do
        sample_name=\$(basename "\$fasta")
        sample_name=\${sample_name%%.*}
        sample_name=\$(echo "\$sample_name" | sed 's/_polished//')
        echo -e "\$(realpath \$fasta)\\t\$sample_name" >> genome_list.txt
    done

    echo "=== Fichier de liste kSNP4 ==="
    cat genome_list.txt
    echo "=== Nombre de génomes : \$(wc -l < genome_list.txt) ==="

    # ─────────────────────────────────────────────────────────────────────────
    # Lancement de kSNP4
    # -in       → fichier de liste
    # -outdir   → dossier de sortie
    # -k        → taille des k-mers (19-21 recommandé pour bactéries)
    # -core     → produit la matrice des SNPs core (présents dans toutes les souches)
    # -min_frac → seuil des SNPs "majoritaires" (présents dans ≥ cette
    #             fraction des souches) ; n'affecte pas les sorties utilisées
    # -CPU      → nombre de threads
    #
    # Le binaire et les options peuvent varier selon la version de kSNP
    # installée : vérifie avec `kSNP4 -h` si tu changes d'image.
    # ─────────────────────────────────────────────────────────────────────────
    kSNP4 \\
        -in genome_list.txt \\
        -outdir ksnp4_results \\
        -k 19 \\
        -core \\
        -min_frac 1.0 \\
        -CPU ${task.cpus}

    # Vérification des sorties et copie à la racine du workdir (requis par
    # les déclarations output:)
    if [ ! -s ksnp4_results/SNPs_all_matrix.fasta ]; then
        echo "ERREUR : kSNP4 n'a pas produit SNPs_all_matrix.fasta" >&2
        echo "Fichiers produits dans ksnp4_results/ :" >&2
        ls -la ksnp4_results/ >&2
        exit 1
    fi
    cp ksnp4_results/SNPs_all_matrix.fasta .

    # Un fichier core peut exister avec des séquences vides (0 SNP core) :
    # on vérifie qu'il contient au moins une ligne de séquence.
    if [ -f ksnp4_results/core_SNPs_matrix.fasta ] && grep -q '^[^>]' ksnp4_results/core_SNPs_matrix.fasta; then
        cp ksnp4_results/core_SNPs_matrix.fasta .
    else
        echo "ATTENTION : aucun SNP core, le rapport core ne sera pas généré" >&2
    fi
    """
}
