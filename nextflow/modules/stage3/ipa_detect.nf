process IPA_DETECT {
    tag "${sample_id}:${chr}"
    publishDir "${params.outdir}/stage3/ipa_detect/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(bam), path(bai), val(chr)
    path intron_annotation
    path rpkm_rds

    output:
    tuple val(sample_id), val(chr), path("cov_${sample_id}_${chr}.RDS"),            emit: coverage
    tuple val(sample_id), val(chr), path("retain_introns_${sample_id}_${chr}.RDS"), emit: retention

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/1_ipa_detect.r \\
        --sample ${sample_id} \\
        --bam ${bam} \\
        --chr ${chr} \\
        --annotation ${intron_annotation} \\
        --rpkm ${rpkm_rds} \\
        --outdir ./
    """
}
