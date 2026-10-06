process PHYLO_REPORT {
    tag "${report_name}"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(report_name), path(fasta), path(treefile)
    path kraken_files
    path ref_kraken_files
    path mlst_files
    path qualimap_dirs
    path checkm2_files

    output:
    path "${report_name}.html"

    script:
    // Nom de la référence = nom du fichier jusqu'au premier point (comme dans main.nf)
    def reference_name = params.reference_fasta ? file(params.reference_fasta).getSimpleName() : ''
    """
    python3 ${projectDir}/bin/phylo_report.py \\
        --fasta        ${fasta} \\
        --tree         ${treefile} \\
        --output       ${report_name}.html \\
        --title        "Comparaison génomique des souches bactériennes" \\
        --kraken_dir   ./ \\
        --mlst_dir     ./ \\
        --qualimap_dir ./ \\
        --checkm2_dir  ./ \\
        --reference_name "${reference_name}" \\
        --pipeline_version "${workflow.manifest.version}"
    """
}
