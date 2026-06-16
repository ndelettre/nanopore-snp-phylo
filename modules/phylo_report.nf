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
    tag "${report_name}"
    label 'process_low'

    publishDir "${params.outdir}",              mode: 'copy'

    input:
    tuple val(report_name), path(fasta), path(treefile)

    output:
    path "${report_name}.html"

    script:
    """
    python3 ${projectDir}/bin/phylo_report.py \\
        --fasta  ${fasta} \\
        --tree   ${treefile} \\
        --output ${report_name}.html \\
        --title  "Phylogénie — ${report_name}"
    """
}
