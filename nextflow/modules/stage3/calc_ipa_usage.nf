process CALC_IPA_USAGE {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage3/calc_ipa_usage/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(bam), path(bai)
    path atlas_rds
    path intron_annotation

    output:
    tuple val(sample_id), path("${params.atlas_name}_${sample_id}_ipa_usage_atlas.csv"), emit: usage

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/8_calculate_ipa_usage_combined.R \\
        --sample ${sample_id} \\
        --bam ${bam} \\
        --atlas ${atlas_rds} \\
        --annotation ${intron_annotation} \\
        --atlas_name ${params.atlas_name} \\
        --outdir ./
    """
}
