# nanopore-snp-phylo

Pipeline de génomique comparative bactérienne pour données Oxford Nanopore MinION — compatible EPI2ME.

## Description

Ce pipeline analyse des souches bactériennes séquencées sur MinION (R10.4.1, kit RBK114) et produit :
- Un rapport HTML interactif avec arbre phylogénétique, matrice de distances SNP, contrôle qualité, identification taxonomique et typage MLST
- Un rapport MultiQC agrégeant les métriques de qualité de toutes les étapes

### Workflow

```
CHOPPER  → NANOSTAT ──────────────────────────────────────────────── MULTIQC
         → KRAKEN2 → BRACKEN (identification taxonomique)                 ↑
         → FLYE → MEDAKA → QUALIMAP ────────────────────────────────────┤
                         → QUAST ────────────────────────────────────────┤
                         → MLST ─────────────────────────────────────────┤
                         → CHECKM2 ──────────────────────────────────────┤
                         → KSNP4 → IQTREE ×2 → PHYLO_REPORT ×2 (tous les SNPs / SNPs core)
```

## Prérequis

- [Nextflow](https://www.nextflow.io/) ≥ 22.10.0
- [Docker](https://www.docker.com/)
- Base de données Kraken2 PlusPF-16 (~15 GB sur disque et en RAM), avec ses fichiers `databaseXXXmers.kmer_distrib` pour Bracken (inclus dans les bases préconstruites)
- Base de données CheckM2 (~2 GB)

### Téléchargement des bases de données

```bash
# Kraken2 PlusPF-16 (version du 26/06/2026 ; nécessite 20 GB de RAM par tâche Kraken2)
mkdir -p /data/kraken2_db
cd /data/kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20260626.tar.gz
tar -xzf k2_pluspf_16_GB_20260626.tar.gz
rm k2_pluspf_16_GB_20260626.tar.gz

# CheckM2
mkdir -p /data/checkm2_db
docker run --rm -v /data/checkm2_db:/data staphb/checkm2:1.0.2 \
    checkm2 database --download --path /data
```

## Installation

### Via EPI2ME Desktop

1. Ouvrir EPI2ME Desktop
2. **Workflows** → **Import Local Workflow** → sélectionner ce dossier
3. Ou importer depuis GitHub : `ndelettre/nanopore-snp-phylo`

### En ligne de commande

```bash
nextflow run ndelettre/nanopore-snp-phylo \
    --fastq_dir /chemin/vers/fastq \
    --kraken_db /data/kraken2_db \
    --checkm2_db /data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
    --outdir output
```

## Paramètres

### Entrées / Sorties

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `fastq_dir` | Dossier contenant les FASTQ (un fichier par barcode) | **obligatoire** |
| `kraken_db` | Dossier de la base Kraken2 | `/data/kraken2_db` |
| `checkm2_db` | Chemin vers le fichier `uniref100.KO.1.dmnd` | `/data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` |
| `outdir` | Dossier de sortie des rapports HTML | `output` |
| `resultsdir` | Dossier de sortie des fichiers intermédiaires | `results` |
| `reference_fasta` | Optionnel : assemblage de référence (FASTA non compressé) ajouté au MLST et à kSNP4. Le nom du fichier jusqu'au premier point devient le nom de la souche | — |

### Filtrage des reads

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `min_length` | Longueur minimale des reads (bp) | `200` |
| `min_quality` | Score de qualité minimum (Phred) | `10` |

### Assemblage et polissage

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `genome_size` | Taille estimée du génome pour Flye | `5m` |
| `medaka_model` | Modèle Medaka selon la flowcell et le basecalling | `r1041_e82_400bps_sup_v5.2.0` |

Modèles Medaka disponibles :

| Flowcell | Kit | Mode | Modèle |
|----------|-----|------|--------|
| R10.4.1 | Kit 14 | SUP | `r1041_e82_400bps_sup_v5.2.0` |
| R10.4.1 | Kit 14 | SUP | `r1041_e82_400bps_sup_v5.0.0` |
| R10.4.1 | Kit 14 | HAC | `r1041_e82_400bps_hac_v5.2.0` |
| R10.4.1 | Kit 14 | SUP | `r1041_e82_400bps_sup_v4.3.0` |
| R10.4.1 | Kit 14 | HAC | `r1041_e82_400bps_hac_v4.3.0` |

### Phylogénie

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `bootstrap` | Activer l'ultrafast bootstrap IQ-TREE (désactivé automatiquement en dessous de 4 souches) | `false` |

## Sorties

```
output/                          ← Rapports visibles dans EPI2ME
├── report_all_snps.html         ← Rapport phylogénétique (tous les SNPs)
├── report_core_snps.html        ← Rapport phylogénétique (SNPs core)
├── multiqc_report.html          ← Rapport QC agrégé
├── pipeline_report.html         ← Métriques d'exécution Nextflow
└── pipeline_timeline.html       ← Chronologie d'exécution

results/                         ← Fichiers intermédiaires
├── chopper/                     ← Reads filtrés
├── nanostat/                    ← Statistiques qualité des reads
├── kraken2/                     ← Rapports taxonomiques
├── bracken/                     ← Abondances réestimées au rang espèce
├── flye/                        ← Assemblages bruts
├── medaka/                      ← Assemblages polishés
├── qualimap/                    ← QC du mapping BAM
├── quast/                       ← QC structurel des assemblages
├── mlst/                        ← Typage MLST
├── checkm2/                     ← Complétude et contamination
├── ksnp4/                       ← Alignements SNP
└── iqtree/                      ← Arbres phylogénétiques (un par rapport)
```

### Rapport phylogénétique

Chaque rapport HTML contient :
- **Contrôle qualité** — tableau par souche avec statut OK/NOK selon :
  - Couverture ≥95% du génome à ≥30X (Qualimap, reads remappés sur l'assemblage poli)
  - Complétude ≥99% (CheckM2)
  - Contamination <1% (CheckM2)
- **Identification taxonomique** — espèces >1% par barcode (Kraken2 + Bracken). Bracken réattribue à l'espèce les reads que Kraken2 laisse au rang genre ; il utilise la plus grande longueur de read disponible dans la base (300 pb pour les bases préconstruites)
- **Typage MLST** — séquence type et allèles (schéma auto-détecté via PubMLST)
- **Arbre phylogénétique** — ML enraciné au midpoint, avec valeurs d'ultrafast bootstrap si `--bootstrap` est activé (fiables à partir de 95)
- **Matrice de distances SNP** — distances pairwise en nombre de SNPs. Les statistiques (min / moyenne / max) portent sur les échantillons du run, sans la référence

Le rapport est autonome : il s'ouvre et s'affiche sans connexion internet.

Deux rapports sont produits, chacun avec son propre arbre IQ-TREE :
- `report_all_snps.html` — basé sur tous les SNPs détectés par kSNP4 (`SNPs_all_matrix.fasta`). Les positions manquantes sont ignorées paire par paire : c'est la référence pour les distances entre souches, notamment avec une souche de référence éloignée
- `report_core_snps.html` — basé sur les SNPs core (présents dans toutes les souches). Plus conservateur, mais une souche éloignée réduit le génome commun et peut sous-estimer les distances entre les autres souches. Il n'est pas produit s'il n'y a aucun SNP core

## Logiciels utilisés

| Outil | Version | Usage |
|-------|---------|-------|
| Chopper | 0.9.0 | Filtrage qualité des reads |
| NanoStat | 1.6.0 | Statistiques des reads |
| Kraken2 | 2.1.3 | Identification taxonomique |
| Bracken | 3.1 | Réestimation des abondances au rang espèce |
| Flye | 2.9.6 | Assemblage de novo |
| Medaka | 2.2.1 | Polissage des assemblages |
| Qualimap | 2.3 | QC du mapping BAM |
| QUAST | 5.2.0 | QC structurel des assemblages |
| MLST | 2.23.0 | Typage MLST |
| CheckM2 | 1.0.2 | Complétude et contamination |
| kSNP4 | 4.0 | SNP calling multi-souches |
| IQ-TREE2 | 2.4.0 | Arbre phylogénétique ML |
| MultiQC | 1.34 | Rapport QC agrégé |

## Notes

- Une souche de référence déjà assemblée peut être ajoutée avec `--reference_fasta`. Elle apparaît dans l'arbre, la matrice SNP et le tableau MLST, avec « N/A » dans le tableau QC. Dans l'identification taxonomique, son espèce est affichée (Kraken2 sur ses contigs, espèce couvrant le plus de paires de bases) avec « N/A » dans les autres colonnes. Elle n'est pas comptée dans le nombre d'échantillons du rapport. Son nom ne doit pas être identique à celui d'un échantillon du run, et elle compte dans les seuils de 3 et 4 souches ci-dessous.
- Le pipeline requiert **au minimum 3 souches** pour produire un arbre phylogénétique, et **au minimum 4** pour le bootstrap (désactivé automatiquement en dessous).
- Le pipeline s'arrête si une souche n'a plus aucun read après filtrage.
- Le modèle Medaka doit correspondre à la flowcell et au mode de basecalling utilisés dans MinKNOW lors du séquençage.
- Les bases de données Kraken2 et CheckM2 sont stockées en dehors du pipeline et réutilisées entre les runs.
- Les assemblages NOK dans le rapport QC ne sont pas exclus automatiquement de la phylogénie — ils sont signalés pour information.

## Auteur

Nicolas Delettre — Hôpitaux Paris Saint-Joseph / Marie-Lannelongue