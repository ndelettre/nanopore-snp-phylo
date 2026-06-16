/*
========================================================================================
    MODULE : kSNP4 — Détection de SNPs et phylogénie sans alignement global
========================================================================================
    kSNP4 identifie les SNPs à partir de k-mers (fragments de séquences de taille k).
    Il ne nécessite PAS de génome de référence et fonctionne sur des génomes complets.

    Sorties principales (dossier -outdir) :
    - SNPs_all                → matrice de TOUS les SNPs (présence/absence par souche)
    - core_SNPs_matrix.fasta  → FASTA des SNPs core (présents dans toutes les souches)
    - tree.core.tre           → arbre NJ basé sur les SNPs core (format Newick)
    - tree.ML.tre             → arbre Maximum Likelihood (format Newick)
    - snp_alignment.fasta     → généré dans ce module à partir de SNPs_all pour IQ-TREE
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
    // pas si ce fichier n'est pas produit. Utile pour les fichiers dont
    // la présence dépend des données (ex: tree.ML.tre si pas assez de SNPs).
    //
    // CORRECTION v2.1 : les noms de fichiers émis ont été corrigés.
    // kSNP4 produit `SNPs_all` (sans extension ni "_matrix"), et les matrices
    // "core" s'appellent `core_SNPs_matrix.fasta`. La cohérence entre les
    // outputs déclarés ici et les `cp` du script est maintenant garantie.
    // ─────────────────────────────────────────────────────────────────────────
    output:
    path "SNPs_all",                                  emit: snp_matrix
    path "core_SNPs_matrix.fasta", optional: true,    emit: core_snp_matrix
    path "snp_alignment.fasta",                       emit: snp_alignment
    path "tree.core.tre",          optional: true,    emit: nj_tree
    path "tree.ML.tre",            optional: true,    emit: ml_tree
    path "ksnp4_results/",                            emit: all_results

    script:
    """
    set -euo pipefail

    # ─────────────────────────────────────────────────────────────────────────
    # kSNP4 requiert un fichier "in_list" au format :
    #   /chemin/absolu/vers/genome.fasta<TAB>NomEchantillon
    # On génère ce fichier dynamiquement à partir des fichiers reçus.
    # basename extrait le nom du fichier, sed supprime le suffixe "_polished".
    # ─────────────────────────────────────────────────────────────────────────
    > genome_list.txt
    for fasta in ${fastas}; do
        sample_name=\$(basename "\$fasta" .fasta | sed 's/_polished//')
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
    # -ML       → calcule un arbre Maximum Likelihood
    # -NJ       → calcule un arbre Neighbor Joining
    # -core     → force le calcul de la matrice/arbre core
    # -min_frac → fraction de souches requise pour qu'un SNP soit "core"
    #             1.0 = 100% (strict), 0.9 = 90% (plus permissif)
    # -CPU      → nombre de threads
    #
    # CORRECTION v2.1 : le binaire et les options peuvent varier selon la
    # version de kSNP installée (kSNP4 vs kSNP4.pl). Le container Biocontainers
    # expose le binaire `kSNP4`. Si tu migres vers une autre version, vérifie
    # avec `kSNP4 -h` les options disponibles.
    # ─────────────────────────────────────────────────────────────────────────
    kSNP4 \\
        -in genome_list.txt \\
        -outdir ksnp4_results \\
        -k 19 \\
        -ML \\
        -NJ \\
        -core \\
        -min_frac 1.0 \\
        -CPU ${task.cpus}

    # ─────────────────────────────────────────────────────────────────────────
    # CORRECTION v2.1 : vérification et copie des sorties selon les VRAIS noms
    # de fichiers produits par kSNP4 (auparavant : "SNPs_all_matrix" qui n'existe
    # pas tel quel — le bon nom est "SNPs_all").
    # ─────────────────────────────────────────────────────────────────────────
    if [ ! -f ksnp4_results/SNPs_all ]; then
        echo "ERREUR : kSNP4 n'a pas produit la matrice SNPs_all"
        echo "Fichiers produits dans ksnp4_results/ :"
        ls -la ksnp4_results/
        exit 1
    fi

    # Copie des fichiers à la racine du workdir (requis par les déclarations output:)
    cp ksnp4_results/SNPs_all                   .
    cp ksnp4_results/core_SNPs_matrix.fasta     . 2>/dev/null || echo "Pas de core_SNPs_matrix.fasta"
    cp ksnp4_results/tree.core.tre              . 2>/dev/null || echo "Pas d'arbre NJ core"
    cp ksnp4_results/tree.ML.tre                . 2>/dev/null || echo "Pas d'arbre ML"

    # ─────────────────────────────────────────────────────────────────────────
    # LEÇON : kSNP4 produit "core_SNPs_matrix.fasta" qui EST déjà un alignement
    # FASTA utilisable par IQ-TREE. On peut donc simplement le copier comme
    # snp_alignment.fasta au lieu de le régénérer.
    #
    # CORRECTION v2.1 : si l'alignement core existe, on l'utilise directement
    # (plus fiable que notre parsing maison). Sinon (cas où -core n'a pas
    # produit de sortie), on retombe sur la génération Python depuis SNPs_all.
    # ─────────────────────────────────────────────────────────────────────────
    if [ -f ksnp4_results/core_SNPs_matrix.fasta ]; then
        echo "Utilisation directe de core_SNPs_matrix.fasta comme alignement"
        cp ksnp4_results/core_SNPs_matrix.fasta snp_alignment.fasta
    else
        echo "Pas de core_SNPs_matrix.fasta, génération depuis SNPs_all"
        python3 << 'PYEOF'
import sys

matrix_file = "SNPs_all"
output_fasta = "snp_alignment.fasta"

try:
    with open(matrix_file) as f:
        lines = f.read().splitlines()

    if not lines:
        print("ERREUR : SNPs_all est vide")
        sys.exit(1)

    # NOTE : Le format de SNPs_all dans kSNP4 est :
    #   Ligne 1        : noms des souches séparés par des espaces/tabs
    #   Lignes suivantes : un nucléotide par souche (colonne) par SNP (ligne)
    # À adapter si ta version de kSNP4 utilise un format différent.
    header = lines[0].strip().split()
    n_samples = len(header)
    sequences = {name: [] for name in header}

    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) == n_samples:
            for i, name in enumerate(header):
                sequences[name].append(parts[i])

    if not any(sequences[n] for n in header):
        print("ERREUR : Aucun SNP dans la matrice")
        sys.exit(1)

    with open(output_fasta, 'w') as out:
        for name in header:
            seq = ''.join(sequences[name])
            out.write(f">{name}\\n{seq}\\n")

    print(f"Alignement FASTA généré : {n_samples} souches, {len(sequences[header[0]])} SNPs")

except FileNotFoundError:
    print(f"ERREUR : Fichier {matrix_file} non trouvé")
    sys.exit(1)
PYEOF
    fi

    # Vérification finale
    if [ ! -s snp_alignment.fasta ]; then
        echo "ERREUR : snp_alignment.fasta n'a pas pu être généré"
        exit 1
    fi
    """
}
