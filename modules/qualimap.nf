process QUALIMAP {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_qualimap/",           emit: results
    path "${sample_id}_qualimapReport.html", emit: html

    script:
    """
    set -euo pipefail
    qualimap bamqc \\
        -bam ${bam} \\
        -outdir ${sample_id}_qualimap \\
        -outformat HTML \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}

    cp ${sample_id}_qualimap/qualimapReport.html ${sample_id}_qualimapReport.html
    """
}
