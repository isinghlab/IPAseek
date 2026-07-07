process COUNT_GENES {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage2/count_genes/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(bam), path(bai)
    path intron_annotation

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_gene_counts.rds"), path("${sample_id}_exon_counts.rds"), emit: counts

    script:
    """
    Rscript ${projectDir}/../../2_gene_preprocessing/1_rnaseq_pipeline/countGeneExpression.R \\
        --sample ${sample_id} \\
        --bam ${bam} \\
        --annotation ${intron_annotation} \\
        --outdir ./ \\
        --paired TRUE
    """
}
