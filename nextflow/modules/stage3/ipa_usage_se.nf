process IPA_USAGE_SE {
    tag 'ipa_usage_se'
    publishDir "${params.outdir}/stage3/ipa_usage_se", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    path usage_files
    val atlas_name

    output:
    path "${atlas_name}_ipa_usage_se.rds", emit: se

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/9_create_ipa_usage_se.R \\
        --atlas_name ${atlas_name} \\
        --outdir ./
    """
}
