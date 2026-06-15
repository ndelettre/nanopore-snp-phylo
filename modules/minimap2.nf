process MINIMAP2 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/minimap2", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    set -euo pipefail
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}
