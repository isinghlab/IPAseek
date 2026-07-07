process FILTER_COVERAGE {
    tag "${sample_id}:${chr}"
    publishDir "${params.outdir}/stage3/filter_coverage/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(chr), path(cov_rds), path(retain_rds), path(retention_summary_rds)

    output:
    tuple val(sample_id), val(chr), path("filtered_${sample_id}_${chr}.RDS"), emit: filtered

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/2_filtering.r \\
        --sample ${sample_id} \\
        --chr ${chr} \\
        --coverage ${cov_rds} \\
        --retention ${retain_rds} \\
        --retention_summary ${retention_summary_rds} \\
        --outdir ./
    """
}
