process PELT_CHANGEPOINT {
    tag "${sample_id}:${chr}"
    publishDir "${params.outdir}/stage3/pelt_changepoint/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(chr), path(filtered_rds), path(cov_rds)

    output:
    tuple val(sample_id), val(chr), path("ipa_final_${sample_id}_filtered_${chr}.RDS"), emit: pelt

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/3_pelt.r \\
        --sample ${sample_id} \\
        --chr ${chr} \\
        --filtered ${filtered_rds} \\
        --coverage ${cov_rds} \\
        --pen_value "c(100,10000)" \\
        --minseglen 200 \\
        --outdir ./
    """
}
