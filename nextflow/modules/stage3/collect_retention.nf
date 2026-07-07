process COLLECT_RETENTION {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage3/collect_retention/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), path(retention_files)

    output:
    tuple val(sample_id), path("retain_all_${sample_id}.RDS"),   emit: retention_all
    tuple val(sample_id), path("retention_data_${sample_id}.rds"), emit: retention_data

    script:
    def ipa_script = "${projectDir}/../../3_ipa_run/scripts/1_ipa_detect.r"
    """
    Rscript - <<'REOF'
    source("${ipa_script}")
    files  <- list.files(".", pattern = "retain_introns_.*\\\\.RDS\$", full.names = TRUE)
    result <- retrieve_intronreten_data(files)
    saveRDS(result[["all"]],            "retain_all_${sample_id}.RDS")
    saveRDS(result[["retention_data"]], "retention_data_${sample_id}.rds")
    REOF
    """
}
