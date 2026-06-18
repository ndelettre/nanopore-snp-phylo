process PHYLO_REPORT {
    tag "${report_name}"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(report_name), path(fasta), path(treefile)
    path kraken_files
    path mlst_files
    path qualimap_dirs
    path checkm2_files

    output:
    path "${report_name}.html"

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
