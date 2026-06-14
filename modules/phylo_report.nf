/*
========================================================================================
    MODULE : PHYLO_REPORT — Rapport HTML interactif (arbre + matrice SNP)
========================================================================================
    Génère deux rapports HTML Plotly :
      - report_all_snps.html  → basé sur snp_alignment.fasta (tous les SNPs)
      - report_core_snps.html → basé sur core_SNPs_matrix.fasta (SNPs core)

    Les deux rapports sont copiés dans output/ pour être détectés par EPI2ME.
========================================================================================
*/
process PHYLO_REPORT {
    label 'process_low'

    publishDir "${params.outdir}/phylo_report", mode: 'copy'
    publishDir "${params.outdir}",              mode: 'copy', pattern: "report_all_snps.html"

    input:
    path treefile
    path snp_alignment
    path core_snps      // peut être un fichier vide si optional

    output:
    path "report_all_snps.html",  emit: report_all
    path "report_core_snps.html", emit: report_core, optional: true

    script:
    """
    set -euo pipefail

    # ── Rapport 1 : tous les SNPs ────────────────────────────────────
    python3 ${projectDir}/bin/phylo_report.py \\
        --fasta  ${snp_alignment} \\
        --tree   ${treefile} \\
        --output report_all_snps.html \\
        --title  "Phylogénie — Tous les SNPs"

    # ── Rapport 2 : SNPs core (seulement si le fichier est non vide) ─
    if [ -s ${core_snps} ]; then
        python3 ${projectDir}/bin/phylo_report.py \\
            --fasta  ${core_snps} \\
            --tree   ${treefile} \\
            --output report_core_snps.html \\
            --title  "Phylogénie — SNPs core"
    else
        echo "Pas de SNPs core disponibles, rapport core ignoré"
    fi
    """
}
