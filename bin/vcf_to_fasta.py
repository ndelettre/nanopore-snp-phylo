#!/usr/bin/env python3
"""
vcf_to_fasta.py
---------------
Convertit un VCF multi-sample (produit par bcftools merge) en alignement
FASTA multi-séquences exploitable par IQ-TREE et phylo_report.py.

Principe :
  - Chaque souche → une séquence FASTA
  - Chaque position SNP variable entre les souches → une colonne de l'alignement
  - Génotypes manquants (./.) → exclusion de la position (core SNPs uniquement)
  - Positions non variables (toutes souches identiques) → exclues
    (pas d'information phylogénétique)
  - Positions FILTER != PASS → exclues
  - Indels → exclus (SNPs uniquement, REF et ALT d'un seul nucléotide)

Approche core SNPs : seules les positions appelées dans TOUTES les souches
sont conservées. Une position avec au moins un génotype manquant (N) est
exclue — garantit des distances pairwise comparables entre toutes les souches.

Le FASTA produit est fonctionnellement équivalent au core_SNPs_matrix.fasta
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
                  variable entre au moins deux souches, appelée dans TOUTES
                  les souches (core SNPs)
    """
    samples = []
    records = []
    n_total = 0
    n_excluded_N = 0
    n_excluded_invariant = 0

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

            ref  = cols[3]            # Allèle de référence
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

            n_total += 1

            # Table des allèles : index 0 = REF, 1 = ALT1, 2 = ALT2, etc.
            alleles = [ref] + alts

            # Index du champ GT (génotype) dans le FORMAT
            gt_idx = fmt.index('GT') if 'GT' in fmt else 0

            # Extraction de la base pour chaque souche
            position = {}
            for i, sample in enumerate(samples):
                sample_data = cols[9 + i].split(':')
                gt_raw = sample_data[gt_idx] if gt_idx < len(sample_data) else '.'

                # Pour les bactéries haploïdes : exclure les génotypes hétérozygotes
                # (artefacts de séquençage sur organismes haploïdes)
                alleles_gt = gt_raw.replace('|', '/').split('/')

                if alleles_gt[0] in ('.', ''):
                    # Génotype manquant
                    position[sample] = 'N'
                elif len(set(alleles_gt)) > 1:
                    # Génotype hétérozygote (0/1, 1/2...) → artefact sur bactérie
                    position[sample] = 'N'
                else:
                    try:
                        allele_idx = int(alleles_gt[0])
                        base = alleles[allele_idx] if allele_idx < len(alleles) else 'N'
                        position[sample] = base
                    except ValueError:
                        position[sample] = 'N'

            # ── Core SNPs : exclure toute position avec au moins un N ────────
            # Garantit que toutes les distances pairwise sont calculées sur
            # exactement les mêmes positions → distances comparables entre
            # toutes les paires de souches.
            if any(b == 'N' for b in position.values()):
                n_excluded_N += 1
                continue

            # Exclure les positions non variables (toutes souches identiques)
            bases = set(position.values())
            if len(bases) == 1:
                n_excluded_invariant += 1
                continue

            records.append(position)

    print(
        f"[vcf_to_fasta] Positions traitées    : {n_total}\n"
        f"[vcf_to_fasta] Exclues (N/manquant)  : {n_excluded_N}\n"
        f"[vcf_to_fasta] Exclues (invariantes) : {n_excluded_invariant}\n"
        f"[vcf_to_fasta] Core SNPs retenus      : {len(records)}"
    )

    return samples, records


# ─────────────────────────────────────────────────────────────────────────────
# 2. ÉCRITURE DU FASTA
# ─────────────────────────────────────────────────────────────────────────────

def write_fasta(samples, records, output_path):
    """
    Écrit l'alignement FASTA multi-séquences.

    Format :
        >souche_A
        ACGT...
        >souche_B
        ACGT...

    Chaque séquence = concaténation des bases aux positions core SNP.
    Lignes de 80 caractères (convention FASTA standard).
    """
    if not records:
        print(
            "[ERREUR] Aucun core SNP trouvé dans le VCF.\n"
            "Vérifiez que le VCF n'est pas vide et que les filtres PASS sont corrects.\n"
            "Si une souche a une couverture insuffisante, elle peut exclure toutes les positions.",
            file=sys.stderr
        )
        sys.exit(1)

    with open(output_path, 'w') as f:
        for sample in samples:
            seq = ''.join(rec[sample] for rec in records)
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
            "pour IQ-TREE et phylo_report.py. Approche core SNPs : seules les\n"
            "positions appelées dans toutes les souches sont conservées."
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
