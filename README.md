# nanopore-snp-phylo

Pipeline SNP et phylogénie pour Nanopore MinION - compatible EPI2ME.

## Description

Analyse de souches bactériennes par séquençage Nanopore :
NanoFilt → NanoStat → Flye → Medaka → kSNP4 → Snippy → IQ-TREE → MultiQC

## Paramètres

- `fastq_dir` : dossier contenant les fichiers FASTQ
- `reference` : génome de référence FASTA (optionnel, requis pour Snippy)
- `outdir` : dossier de sortie (défaut: results)
- `min_length` : longueur minimale des reads (défaut: 200)
- `min_quality` : qualité minimale Phred (défaut: 10)
- `genome_size` : taille estimée du génome pour Flye (défaut: 5m)
- `medaka_model` : modèle Medaka selon la flowcell (défaut: r941_min_high_g360)
