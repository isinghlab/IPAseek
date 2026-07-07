process INTRON_RETENTION_SE {
    tag 'intron_retention_se'
    publishDir "${params.outdir}/stage3/intron_retention_se", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    path retention_files

    output:
    path "se_intron_retention.RDS", emit: se

    script:
    def ipa_script = "${projectDir}/../../3_ipa_run/scripts/1_ipa_detect.r"
    """
    Rscript - <<'REOF'
    source("${ipa_script}")
    files <- list.files(".", pattern = "retention_data_.*\\\\.rds\$", full.names = TRUE)
    se    <- intronret_se(files)
    saveRDS(se, "se_intron_retention.RDS")
    REOF
    """
}
