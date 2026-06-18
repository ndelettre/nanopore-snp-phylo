/*
========================================================================================
    MODULE : PHYLO_REPORT — Rapport HTML interactif
========================================================================================
    Génère un rapport HTML complet par analyse SNP :
      - Header HPSJ + titre + date
      - Contrôle qualité des assemblages (Qualimap + CheckM2)
      - Identification taxonomique (Kraken2)
      - Typage MLST
      - Arbre phylogénétique IQ-TREE (Plotly)
      - Matrice de distances SNP (Plotly)
      - Footer avec versions des logiciels
========================================================================================
*/
process PHYLO_REPORT {
    tag "${report_name}"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
tuple val(report_name), path(fasta), path(treefile)
path kraken_files    // liste des *.kraken2.report
path mlst_files      // liste des *.mlst.tsv
path qualimap_dirs   // liste des *_qualimap/
path checkm2_files   // liste des *_checkm2_quality_report.tsv

script:
"""
python3 ${projectDir}/bin/phylo_report.py \\
    --fasta        ${fasta} \\
    --tree         ${treefile} \\
    --output       ${report_name}.html \\
    --title        "Comparaison génomique des souches bactériennes" \\
    --kraken_dir   ./ \\
    --mlst_dir     ./ \\
    --qualimap_dir ./ \\
    --checkm2_dir  ./
"""

}
