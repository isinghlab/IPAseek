process FILTER_CHANGEPOINTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage3/filter_changepoints/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(pelt_files), path(bam), path(bai)
    path intron_annotation
    path rpkm_rds

    output:
    tuple val(sample_id), path("filter_cpts_tpm_exprs*_${sample_id}.csv"), path("filter_cpts_all_${sample_id}.csv"), emit: cpts

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/4_filter_cpts_de.R \\
        --sample ${sample_id} \\
        --bam ${bam} \\
        --annotation ${intron_annotation} \\
        --rpkm ${rpkm_rds} \\
        --outdir ./
    """
}
