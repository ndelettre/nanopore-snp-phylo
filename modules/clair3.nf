process CLAIR3 {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/clair3/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail
    samtools faidx ${reference}

    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=/opt/models/r1041_e82_400bps_sup_v420 \\
        --output=clair3_out \\
        --include_all_ctgs \\
        --no_phasing_for_fa

    cp clair3_out/merge_output.vcf.gz ${sample_id}.vcf.gz
    """
}
