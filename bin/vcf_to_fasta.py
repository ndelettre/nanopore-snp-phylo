#!/usr/bin/env python3
"""
vcf_to_fasta.py
---------------
Convertit un VCF multi-sample (bcftools merge) en alignement FASTA multi-séquences.
Chaque séquence = une souche, chaque position = un SNP variable entre les souches.
Les positions avec données manquantes (./.) sont remplacées par 'N'.

Usage :
    python3 vcf_to_fasta.py \
        --vcf  merged.vcf.gz \
        --output snp_alignment.fasta

Dépendances : pysam
"""

import argparse
import gzip
import sys


def parse_vcf(vcf_path):
    """
    Lit un VCF multi-sample gzippé ou non.
    Retourne :
        samples : liste des noms de souches (dans l'ordre des colonnes)
        records : liste de dict {sample: base} pour chaque position SNP
    """
    samples = []
    records = []

    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # Ligne d'en-tête : colonnes 9+ = noms des souches
                cols = line.strip().split('\t')
                samples = cols[9:]
                continue

            cols = line.strip().split('\t')
            if len(cols) < 9 + len(samples):
                continue

            ref  = cols[3]
            alts = cols[4].split(',')
            fmt  = cols[8].split(':')

            # On ne garde que les SNPs bi-alléliques (REF et ALT d'un seul nucléotide)
            if len(ref) != 1:
                continue
            if not all(len(a) == 1 for a in alts):
                continue
            # Exclure les positions FILTER != PASS (sauf '.')
            filt = cols[6]
            if filt not in ('PASS', '.'):
                continue

            alleles = [ref] + alts  # index 0 = REF, 1 = ALT1, etc.

            gt_idx = fmt.index('GT') if 'GT' in fmt else 0

            position = {}
            is_variable = False

            for i, sample in enumerate(samples):
                sample_data = cols[9 + i].split(':')
                gt_raw = sample_data[gt_idx] if gt_idx < len(sample_data) else '.'

                # Normaliser GT : 0/0, 0|0, 1/1, ./. etc.
                gt = gt_raw.replace('|', '/').split('/')[0]

                if gt == '.' or gt == '':
                    position[sample] = 'N'
                else:
                    try:
                        allele_idx = int(gt)
                        base = alleles[allele_idx] if allele_idx < len(alleles) else 'N'
                        position[sample] = base
                    except ValueError:
                        position[sample] = 'N'

            # Vérifier que la position est réellement variable entre les souches
            bases = set(b for b in position.values() if b != 'N')
            if len(bases) > 1:
                is_variable = True

            if is_variable:
                records.append(position)

    return samples, records


def write_fasta(samples, records, output_path):
    """
    Écrit l'alignement FASTA multi-séquences.
    Une séquence par souche = concaténation de toutes les bases aux positions SNP.
    """
    if not records:
        print("[ERREUR] Aucun SNP variable trouvé dans le VCF.", file=sys.stderr)
        sys.exit(1)

    print(f"[vcf_to_fasta] {len(samples)} souches, {len(records)} positions SNP variables")

    with open(output_path, 'w') as f:
        for sample in samples:
            seq = ''.join(rec.get(sample, 'N') for rec in records)
            f.write(f'>{sample}\n')
            # Écrire en lignes de 80 caractères
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="VCF multi-sample → alignement FASTA SNP pour IQ-TREE / phylo_report"
    )
    parser.add_argument('--vcf',    required=True, help='VCF multi-sample (gzippé ou non)')
    parser.add_argument('--output', required=True, help='Fichier FASTA de sortie')
    args = parser.parse_args()

    print(f"[vcf_to_fasta] Lecture : {args.vcf}")
    samples, records = parse_vcf(args.vcf)

    if not samples:
        print("[ERREUR] Aucun sample trouvé dans le VCF.", file=sys.stderr)
        sys.exit(1)

    write_fasta(samples, records, args.output)
    print(f"[vcf_to_fasta] FASTA écrit : {args.output}")


if __name__ == '__main__':
    main()
