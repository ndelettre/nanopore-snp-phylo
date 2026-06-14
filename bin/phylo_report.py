#!/usr/bin/env python3
"""
phylo_report.py
---------------
Génère un rapport HTML interactif contenant :
  - Un arbre phylogénétique IQ-TREE enraciné au midpoint
  - Une matrice de distances SNP pairwise calculée depuis un alignement FASTA

Les deux visualisations sont ordonnées de manière identique (ordre visuel
de l'arbre, de haut en bas) pour faciliter la lecture.

Usage :
    python3 phylo_report.py \
        --fasta  snp_alignment.fasta \
        --tree   phylo_tree.treefile \
        --output report.html \
        --title  "Rapport SNP — HPSJ"

Dépendances : biopython, plotly, pandas, numpy
"""

import argparse
import sys
from itertools import combinations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from copy import deepcopy


# ─────────────────────────────────────────────────────────────────────────────
# 1. LECTURE DU FASTA ET CALCUL DES DISTANCES SNP
# ─────────────────────────────────────────────────────────────────────────────

def read_fasta(path):
    """
    Lit un alignement FASTA.
    Retourne un dict {nom_sequence: str_sequence}.
    """
    sequences = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]   # prendre uniquement le premier mot
                sequences[current] = []
            elif current:
                sequences[current].append(line)
    return {k: "".join(v) for k, v in sequences.items()}


def compute_snp_distances(sequences):
    """
    Calcule la matrice de distances SNP pairwise.
    Distance = nombre de positions différentes entre deux séquences alignées.
    Ignore les positions avec N ou - dans l'une ou l'autre séquence.

    Retourne un DataFrame carré symétrique (diagonale = 0).
    """
    samples = list(sequences.keys())
    mat = pd.DataFrame(0, index=samples, columns=samples, dtype=int)

    for s1, s2 in combinations(samples, 2):
        seq1 = sequences[s1].upper()
        seq2 = sequences[s2].upper()
        dist = sum(
            a != b
            for a, b in zip(seq1, seq2)
            if a not in "N-" and b not in "N-"
        )
        mat.loc[s1, s2] = dist
        mat.loc[s2, s1] = dist

    return mat


# ─────────────────────────────────────────────────────────────────────────────
# 2. CHARGEMENT ET ENRACINEMENT DE L'ARBRE
# ─────────────────────────────────────────────────────────────────────────────

def load_and_root_tree(tree_path):
    """
    Charge un arbre Newick IQ-TREE et l'enracine au midpoint.

    L'enracinement au midpoint place la racine au milieu de la plus longue
    branche de l'arbre — méthode standard quand on n'a pas d'outgroup.

    Retourne l'arbre enraciné (Bio.Phylo.BaseTree.Tree).
    """
    tree = Phylo.read(tree_path, "newick")
    tree.root_at_midpoint()
    return tree


def has_bootstrap(tree):
    """Détecte si l'arbre contient des valeurs de bootstrap."""
    return any(
        c.confidence is not None
        for c in tree.find_clades()
        if not c.is_terminal()
    )


# ─────────────────────────────────────────────────────────────────────────────
# 3. CONVERSION ARBRE → COORDONNÉES PLOTLY
# ─────────────────────────────────────────────────────────────────────────────

def get_tree_coords(tree):
    """
    Calcule les coordonnées (x, y) de chaque nœud pour le rendu Plotly.

    Convention :
      x = distance cumulée depuis la racine (longueur de branche)
      y = position verticale :
            feuilles  → espacées régulièrement 0, 1, 2, ...
            internes  → centrés entre leurs enfants (haut + bas) / 2

    Retourne :
        nodes  : dict {clade → (x, y)}
        leaves : list de clades terminaux dans l'ordre visuel (haut → bas)
    """
    leaves = tree.get_terminals()
    leaf_y = {leaf: i for i, leaf in enumerate(leaves)}

    nodes = {}

    def assign_coords(clade, x_start):
        branch_length = clade.branch_length if clade.branch_length else 0.0
        x = x_start + branch_length

        if clade.is_terminal():
            y = leaf_y[clade]
        else:
            for child in clade.clades:
                assign_coords(child, x)
            child_ys = [nodes[child][1] for child in clade.clades]
            y = (min(child_ys) + max(child_ys)) / 2

        nodes[clade] = (x, y)

    tree.root.branch_length = 0.0
    assign_coords(tree.root, 0.0)

    return nodes, leaves


def build_tree_traces(tree, nodes, bootstrap=True):
    """
    Construit les traces Plotly pour l'arbre :
      - Lignes de branches (horizontales + verticales)
      - Marqueurs sur les feuilles avec noms
      - Carrés bootstrap colorés sur les nœuds internes
    """
    traces = []

    # ── Branches ────────────────────────────────────────────────────
    edge_x, edge_y = [], []

    def add_edges(clade):
        px, py = nodes[clade]
        for child in clade.clades:
            cx, cy = nodes[child]
            edge_x.extend([px, px, None])   # ligne verticale
            edge_y.extend([cy, py, None])
            edge_x.extend([px, cx, None])   # ligne horizontale
            edge_y.extend([cy, cy, None])
            add_edges(child)

    add_edges(tree.root)

    traces.append(go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line=dict(color="#2c3e50", width=1.8),
        hoverinfo="skip",
        showlegend=False
    ))

    # ── Feuilles ────────────────────────────────────────────────────
    leaf_x     = [nodes[c][0] for c in tree.get_terminals()]
    leaf_y_val = [nodes[c][1] for c in tree.get_terminals()]
    leaf_names = [c.name      for c in tree.get_terminals()]

    traces.append(go.Scatter(
        x=leaf_x, y=leaf_y_val,
        mode="markers+text",
        text=leaf_names,
        textposition="middle right",
        textfont=dict(size=11, color="#2c3e50"),
        marker=dict(size=7, color="#e74c3c"),
        hovertemplate="%{text}<extra></extra>",
        showlegend=False
    ))

    # ── Bootstrap ───────────────────────────────────────────────────
    if bootstrap:
        bs_x, bs_y, bs_text, bs_colors = [], [], [], []
        for clade in tree.find_clades(order="level"):
            if not clade.is_terminal() and clade.confidence is not None:
                x, y = nodes[clade]
                val  = clade.confidence
                bs_x.append(x)
                bs_y.append(y)
                bs_text.append(str(int(val)))
                bs_colors.append(
                    "#27ae60" if val >= 90 else
                    "#f39c12" if val >= 70 else
                    "#e74c3c"
                )

        if bs_x:
            traces.append(go.Scatter(
                x=bs_x, y=bs_y,
                mode="markers+text",
                text=bs_text,
                textposition="top center",
                textfont=dict(size=9),
                marker=dict(size=6, color=bs_colors, symbol="square"),
                hovertemplate="Bootstrap : %{text}<extra></extra>",
                showlegend=False
            ))

    return traces


# ─────────────────────────────────────────────────────────────────────────────
# 4. RÉORDONNANCEMENT DE LA MATRICE
# ─────────────────────────────────────────────────────────────────────────────

def reorder_matrix(matrix, leaves):
    """
    Réordonne la matrice dans le même ordre visuel que l'arbre (haut → bas).
    Signale les samples présents dans un seul des deux fichiers.
    """
    leaf_names = [c.name for c in leaves]

    missing_from_matrix = [n for n in leaf_names if n not in matrix.index]
    missing_from_tree   = [n for n in matrix.index if n not in leaf_names]

    if missing_from_matrix:
        print(f"[WARN] Dans l'arbre mais absent du FASTA : {missing_from_matrix}",
              file=sys.stderr)
    if missing_from_tree:
        print(f"[WARN] Dans le FASTA mais absent de l'arbre : {missing_from_tree}",
              file=sys.stderr)

    common = [n for n in leaf_names if n in matrix.index]
    return matrix.loc[common, common]


# ─────────────────────────────────────────────────────────────────────────────
# 5. CONSTRUCTION DE LA HEATMAP
# ─────────────────────────────────────────────────────────────────────────────

def build_heatmap_trace(matrix):
    """
    Heatmap Plotly de la matrice de distances SNP.
    Palette : blanc (0 SNP) → bleu foncé (distance max).
    Annotations avec valeurs numériques dans chaque cellule.
    """
    samples = list(matrix.index)
    z       = matrix.values.tolist()

    # Annotations : valeur entière dans chaque cellule.
    # On normalise chaque valeur entre 0 et 1 par rapport au max.
    # Colorscale Blues inversée : blanc (0) → bleu foncé (max).
    # Texte blanc si case foncée (normalisé > 0.6), gris foncé sinon.
    max_val = float(matrix.values.max()) or 1.0  # éviter division par zéro

    annotations = []
    for i, row in enumerate(samples):
        for j, col in enumerate(samples):
            val      = int(matrix.iloc[i, j])
            norm_val = val / max_val
            annotations.append(dict(
                x=col, y=row,
                text=str(val),
                showarrow=False,
                font=dict(
                    size=10,
                    color="white" if norm_val > 0.6 else "#2c3e50"
                )
            ))

    trace = go.Heatmap(
        z=z,
        x=samples,
        y=samples,
        colorscale="Blues",
        reversescale=True,
        colorbar=dict(
            title=dict(text="SNPs", side="right"),
            len=0.6
        ),
        hovertemplate="%{y} vs %{x}<br><b>%{z} SNPs</b><extra></extra>",
        showscale=True
    )

    return trace, annotations


# ─────────────────────────────────────────────────────────────────────────────
# 6. ASSEMBLAGE DU RAPPORT HTML
# ─────────────────────────────────────────────────────────────────────────────

def build_report(fasta_path, tree_path, output_path, title="Rapport Phylogénétique"):

    print(f"[1/6] Lecture du FASTA : {fasta_path}")
    sequences = read_fasta(fasta_path)
    n_samples = len(sequences)
    print(f"      {n_samples} séquences lues")

    print(f"[2/6] Calcul des distances SNP pairwise...")
    matrix_raw = compute_snp_distances(sequences)

    print(f"[3/6] Chargement et enracinement au midpoint : {tree_path}")
    tree = load_and_root_tree(tree_path)
    bootstrap_present = has_bootstrap(tree)
    print(f"      Bootstrap détecté : {bootstrap_present}")

    print(f"[4/6] Calcul des coordonnées de l'arbre...")
    nodes, leaves = get_tree_coords(tree)

    print(f"[5/6] Réordonnancement de la matrice selon l'arbre...")
    matrix = reorder_matrix(matrix_raw, leaves)

    # ── Statistiques ────────────────────────────────────────────────
    dist_vals = matrix.values[np.triu_indices(len(matrix), k=1)]
    dist_min  = int(dist_vals.min())  if len(dist_vals) else 0
    dist_max  = int(dist_vals.max())  if len(dist_vals) else 0
    dist_mean = f"{dist_vals.mean():.1f}" if len(dist_vals) else "N/A"

    # Paires extrêmes
    flat = matrix.where(np.triu(np.ones(matrix.shape, dtype=bool), k=1))
    pair_min = flat.stack().idxmin()
    pair_max = flat.stack().idxmax()

    # ── Figures Plotly ───────────────────────────────────────────────
    print(f"[6/6] Génération du HTML : {output_path}")

    n_leaves       = len(leaves)
    tree_height    = max(500, n_leaves * 38)
    heatmap_size   = max(420, n_leaves * 44)
    max_name_len   = max(len(c.name) for c in leaves)
    right_margin   = max(130, max_name_len * 7)

    # Arbre
    tree_traces = build_tree_traces(tree, nodes, bootstrap=bootstrap_present)
    fig_tree = go.Figure(data=tree_traces)
    fig_tree.update_layout(
        title=dict(text="Arbre Phylogénétique — enraciné au midpoint", x=0.5,
                   font=dict(size=15)),
        height=tree_height,
        xaxis=dict(title="Distance évolutive", showgrid=True,
                   gridcolor="#ecf0f1", zeroline=False),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-0.5, n_leaves - 0.5]),
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(l=20, r=right_margin, t=50, b=50),
        hovermode="closest"
    )

    if bootstrap_present:
        fig_tree.add_annotation(
            text="Bootstrap : <span style='color:#27ae60'>■</span> ≥90  "
                 "<span style='color:#f39c12'>■</span> ≥70  "
                 "<span style='color:#e74c3c'>■</span> <70",
            xref="paper", yref="paper", x=0.01, y=0.01,
            showarrow=False, font=dict(size=10),
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="#bdc3c7", borderwidth=1
        )

    # Heatmap
    heatmap_trace, annotations = build_heatmap_trace(matrix)
    fig_heatmap = go.Figure(data=[heatmap_trace])
    fig_heatmap.update_layout(
        title=dict(text="Matrice de Distances SNP", x=0.5, font=dict(size=15)),
        height=heatmap_size,
        width=heatmap_size + 120,
        annotations=annotations,
        xaxis=dict(tickangle=45, tickfont=dict(size=10), side="bottom"),
        yaxis=dict(tickfont=dict(size=10), autorange="reversed"),
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(l=110, r=80, t=60, b=130)
    )

    # Export HTML
    tree_html = fig_tree.to_html(
        full_html=False, include_plotlyjs="cdn", div_id="tree_div"
    )
    heatmap_html = fig_heatmap.to_html(
        full_html=False, include_plotlyjs=False, div_id="heatmap_div"
    )

    # ── Paires extrêmes pour le résumé ──────────────────────────────
    pair_min_str = f"{pair_min[0]} / {pair_min[1]} ({dist_min} SNPs)"
    pair_max_str = f"{pair_max[0]} / {pair_max[1]} ({dist_max} SNPs)"

    html = f"""<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: #f5f6fa;
            color: #2c3e50;
        }}
        header {{
            background: linear-gradient(135deg, #1a252f 0%, #2c3e50 100%);
            color: white;
            padding: 22px 40px;
            border-bottom: 4px solid #3498db;
        }}
        header h1 {{ font-size: 1.6em; font-weight: 600; }}
        header p  {{ margin-top: 5px; opacity: 0.65; font-size: 0.85em; }}

        .container {{ max-width: 1400px; margin: 0 auto; padding: 24px 32px; }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 14px;
            margin-bottom: 26px;
        }}
        .stat-card {{
            background: white;
            border-radius: 10px;
            padding: 16px 18px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.07);
            border-left: 4px solid #3498db;
        }}
        .stat-card.green  {{ border-left-color: #27ae60; }}
        .stat-card.orange {{ border-left-color: #e67e22; }}
        .stat-card.red    {{ border-left-color: #e74c3c; }}
        .stat-card .value {{
            font-size: 1.7em; font-weight: 700; color: #2980b9;
        }}
        .stat-card.green  .value {{ color: #27ae60; }}
        .stat-card.orange .value {{ color: #e67e22; }}
        .stat-card.red    .value {{ color: #e74c3c; }}
        .stat-card .label {{
            font-size: 0.75em; color: #7f8c8d; margin-top: 4px;
            text-transform: uppercase; letter-spacing: 0.05em;
        }}
        .stat-card .sub {{
            font-size: 0.78em; color: #95a5a6; margin-top: 3px;
        }}

        .section {{
            background: white;
            border-radius: 12px;
            padding: 22px;
            margin-bottom: 26px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.07);
        }}
        .section h2 {{
            font-size: 0.85em; color: #95a5a6;
            text-transform: uppercase; letter-spacing: 0.08em;
            margin-bottom: 14px; padding-bottom: 10px;
            border-bottom: 2px solid #ecf0f1;
        }}
        footer {{
            text-align: center; padding: 18px;
            color: #bdc3c7; font-size: 0.78em;
        }}
    </style>
</head>
<body>

<header>
    <h1>{title}</h1>
    <p>Arbre IQ-TREE enraciné au midpoint · Distances calculées depuis l'alignement SNP</p>
</header>

<div class="container">

    <div class="stats-grid">
        <div class="stat-card">
            <div class="value">{n_samples}</div>
            <div class="label">Échantillons</div>
        </div>
        <div class="stat-card green">
            <div class="value">{dist_min}</div>
            <div class="label">Distance minimale</div>
            <div class="sub">{pair_min_str}</div>
        </div>
        <div class="stat-card orange">
            <div class="value">{dist_mean}</div>
            <div class="label">Distance moyenne (SNPs)</div>
        </div>
        <div class="stat-card red">
            <div class="value">{dist_max}</div>
            <div class="label">Distance maximale</div>
            <div class="sub">{pair_max_str}</div>
        </div>
        <div class="stat-card {'green' if bootstrap_present else 'red'}">
            <div class="value" style="font-size:1.1em;">
                {'✓ Présent' if bootstrap_present else '✗ Absent'}
            </div>
            <div class="label">Bootstrap IQ-TREE</div>
        </div>
    </div>

    <div class="section">
        <h2>Arbre Phylogénétique</h2>
        {tree_html}
    </div>

    <div class="section">
        <h2>Matrice de Distances SNP</h2>
        {heatmap_html}
    </div>

</div>

<footer>
    Généré avec phylo_report.py · Biopython · Plotly
</footer>

</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"\n✓ Rapport généré : {output_path}")
    print(f"  {n_samples} samples | SNPs : {dist_min}–{dist_max} (moy. {dist_mean})")
    print(f"  Paire la plus proche : {pair_min_str}")
    print(f"  Paire la plus distante : {pair_max_str}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Rapport HTML phylogénétique interactif (arbre midpoint + matrice SNP)"
    )
    parser.add_argument("--fasta",   required=True, help="Alignement SNP FASTA (kSNP4)")
    parser.add_argument("--tree",    required=True, help="Arbre IQ-TREE Newick (.treefile)")
    parser.add_argument("--output",  default="report.html", help="Fichier HTML de sortie")
    parser.add_argument("--title",   default="Rapport Phylogénétique", help="Titre du rapport")
    args = parser.parse_args()

    build_report(
        fasta_path  = args.fasta,
        tree_path   = args.tree,
        output_path = args.output,
        title       = args.title
    )


if __name__ == "__main__":
    main()
