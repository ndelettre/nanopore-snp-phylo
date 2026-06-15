process VCF_TO_FASTA {
    label 'process_low'
    publishDir "${params.outdir}/snp_alignment", mode: 'copy'

    input:
    path merged_vcf

    output:
    path "snp_alignment.fasta", emit: fasta

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/vcf_to_fasta.py \\
        --vcf    ${merged_vcf} \\
        --output snp_alignment.fasta
    """
}
