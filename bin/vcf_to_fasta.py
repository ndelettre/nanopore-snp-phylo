#!/usr/bin/env python3
"""
vcf_to_fasta.py
---------------
Convertit un VCF multi-sample (produit par bcftools merge) en alignement
FASTA multi-séquences exploitable par IQ-TREE et phylo_report.py.

Principe :
  - Chaque souche → une séquence FASTA
  - Chaque position SNP variable entre les souches → une colonne de l'alignement
  - Génotypes manquants (./.) → 'N'
  - Positions non variables (toutes souches identiques) → exclues
    (pas d'information phylogénétique)
  - Positions FILTER != PASS → exclues
  - Indels → exclus (SNPs uniquement, REF et ALT d'un seul nucléotide)

Le FASTA produit est fonctionnellement équivalent au snp_alignment.fasta
de kSNP4 — PHYLO_REPORT et IQ-TREE le consomment de manière identique.

Usage :
    python3 vcf_to_fasta.py \\
        --vcf  merged.vcf.gz \\
        --output snp_alignment.fasta

Dépendances : stdlib Python uniquement (gzip, argparse, sys)
"""

import argparse
import gzip
import sys


# ─────────────────────────────────────────────────────────────────────────────
# 1. LECTURE DU VCF
# ─────────────────────────────────────────────────────────────────────────────

def parse_vcf(vcf_path):
    """
    Lit un VCF multi-sample (gzippé ou non).

    Retourne :
        samples : liste des noms de souches dans l'ordre des colonnes VCF
        records : liste de dict {sample_id: base} pour chaque position SNP
                  variable entre au moins deux souches
    """
    samples = []
    records = []

    # Supporte les VCF compressés (.vcf.gz) et non compressés (.vcf)
    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:

            # Lignes de métadonnées : ignorer
            if line.startswith('##'):
                continue

            # Ligne d'en-tête : extraire les noms des souches (colonnes 9+)
            if line.startswith('#CHROM'):
                cols = line.strip().split('\t')
                samples = cols[9:]
                continue

            # Ligne de variant
            cols = line.strip().split('\t')
            if len(cols) < 9 + len(samples):
                continue

            ref  = cols[3]           # Allèle de référence
            alts = cols[4].split(',') # Allèles alternatifs (peut être multi-allélique)
            fmt  = cols[8].split(':') # Format des champs génotype
            filt = cols[6]            # Statut FILTER

            # Exclusion des indels : on ne garde que les SNPs bi-nucléotidiques
            # (REF = 1 base, tous les ALT = 1 base)
            if len(ref) != 1:
                continue
            if not all(len(a) == 1 for a in alts):
                continue

            # Exclusion des positions non validées
            # '.' = pas de filtre appliqué (accepté), 'PASS' = validé
            if filt not in ('PASS', '.'):
                continue

            # Table des allèles : index 0 = REF, 1 = ALT1, 2 = ALT2, etc.
            alleles = [ref] + alts

            # Index du champ GT (génotype) dans le FORMAT
            gt_idx = fmt.index('GT') if 'GT' in fmt else 0

            # Extraction de la base pour chaque souche
            position = {}
            for i, sample in enumerate(samples):
                sample_data = cols[9 + i].split(':')
                gt_raw = sample_data[gt_idx] if gt_idx < len(sample_data) else '.'

                # Normalisation du génotype : 0/0, 0|0, 1/1, ./. → premier allèle
                # Pour les bactéries haploïdes, le premier allèle suffit
                gt = gt_raw.replace('|', '/').split('/')[0]

                if gt in ('.', ''):
                    # Génotype manquant → N (position non appelée pour cette souche)
                    position[sample] = 'N'
                else:
                    try:
                        allele_idx = int(gt)
                        base = alleles[allele_idx] if allele_idx < len(alleles) else 'N'
                        position[sample] = base
                    except ValueError:
                        position[sample] = 'N'

            # Ne garder que les positions réellement variables entre les souches
            # (au moins deux bases différentes en excluant les N)
            bases = set(b for b in position.values() if b != 'N')
            if len(bases) > 1:
                records.append(position)

    return samples, records


# ─────────────────────────────────────────────────────────────────────────────
# 2. ÉCRITURE DU FASTA
# ─────────────────────────────────────────────────────────────────────────────

def write_fasta(samples, records, output_path):
    """
    Écrit l'alignement FASTA multi-séquences.

    Format :
        >souche_A
        ACGTNN...
        >souche_B
        ACGTAN...

    Chaque séquence = concaténation des bases aux positions SNP variables.
    Lignes de 80 caractères (convention FASTA standard).
    """
    if not records:
        print(
            "[ERREUR] Aucun SNP variable trouvé dans le VCF.\n"
            "Vérifiez que le VCF n'est pas vide et que les filtres PASS sont corrects.",
            file=sys.stderr
        )
        sys.exit(1)

    n_snps = len(records)
    print(f"[vcf_to_fasta] {len(samples)} souches | {n_snps} positions SNP variables")

    with open(output_path, 'w') as f:
        for sample in samples:
            seq = ''.join(rec.get(sample, 'N') for rec in records)
            f.write(f'>{sample}\n')
            # Écriture en lignes de 80 caractères (convention FASTA)
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    print(f"[vcf_to_fasta] FASTA écrit : {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. POINT D'ENTRÉE
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convertit un VCF multi-sample (bcftools merge) en alignement FASTA\n"
            "pour IQ-TREE et phylo_report.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--vcf',
        required=True,
        help='VCF multi-sample produit par bcftools merge (gzippé ou non)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Fichier FASTA de sortie (ex: snp_alignment.fasta)'
    )
    args = parser.parse_args()

    print(f"[vcf_to_fasta] Lecture : {args.vcf}")
    samples, records = parse_vcf(args.vcf)

    if not samples:
        print("[ERREUR] Aucun sample trouvé dans le VCF.", file=sys.stderr)
        sys.exit(1)

    write_fasta(samples, records, args.output)


if __name__ == '__main__':
    main()
