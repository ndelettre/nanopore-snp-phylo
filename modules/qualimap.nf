// modules/qualimap.nf
process QUALIMAP {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qualimap/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_qualimap/",                              emit: results
    path "${sample_id}_qualimap/genome_results.txt",            emit: stats

    script:
    """
    set -euo pipefail
    qualimap bamqc \\
        -bam ${bam} \\
        -outdir ${sample_id}_qualimap \\
        -outformat HTML \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}
