process FILTER_TERMINAL_EXONS {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage3/filter_terminal_exons/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), path(cpt_exprs_all), path(bam), path(bai)
    path intron_annotation

    output:
    tuple val(sample_id), path("te_expression${sample_id}.csv"),            emit: te
    tuple val(sample_id), path("te_expression${sample_id}_unfiltered.csv"), emit: te_unfiltered

    script:
    """
    Rscript ${projectDir}/../../3_ipa_run/scripts/6_analyse_exon_structure.r \\
        --sample ${sample_id} \\
        --cpt_exprs ${cpt_exprs_all} \\
        --bam ${bam} \\
        --annotation ${intron_annotation} \\
        --outdir ./
    """
}
