process MAKE_ATLAS {
    tag 'make_atlas'
    publishDir "${params.outdir}/stage3/make_atlas", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    path te_files
    path intron_annotation
    val atlas_name

    output:
    path "${atlas_name}_full_ipa_atlas_conf.RDS", emit: atlas_rds
    path "${atlas_name}_full_ipa_atlas_conf.csv", emit: atlas_csv

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/7_make_atlas.r \\
        --atlas_name ${atlas_name} \\
        --annotation ${intron_annotation} \\
        --outdir ./
    """
}
