process BCFTOOLS_MERGE {
    label 'process_low'
    publishDir "${params.outdir}/bcftools", mode: 'copy'

    input:
    path vcf_list

    output:
    path "merged.vcf.gz", emit: merged_vcf

    script:
    """
    set -euo pipefail

    for vcf in ${vcf_list}; do
        bcftools index --tbi \$vcf
    done

    bcftools merge \\
        --merge snps \\
        --output-type z \\
        --output merged.vcf.gz \\
        ${vcf_list}
    """
}
