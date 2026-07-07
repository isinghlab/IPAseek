process MERGE_CHANGEPOINTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage3/merge_changepoints/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), path(tpm_exprs_files), path(cpts_all_file)

    output:
    tuple val(sample_id), path("${sample_id}_cpt_all.csv"),       emit: cpt_all
    tuple val(sample_id), path("${sample_id}_cpt_exprs_all.csv"), emit: merged

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/5_merge.r \\
        --sample ${sample_id} \\
        --outdir ./
    """
}
