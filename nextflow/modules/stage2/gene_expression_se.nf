process GENE_EXPRESSION_SE {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage2/gene_expression_se/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(gene_counts), path(exon_counts)
    path intron_annotation

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_rpkm.rds"), emit: rpkm
    tuple val(sample_id), val(condition), path("${sample_id}_all.rds"),  emit: all_se

    script:
    def gene_expr_script = "${projectDir}/../../2_gene_preprocessing/1_rnaseq_pipeline/gene_expr.R"
    """
    Rscript - <<'REOF'
    source("${gene_expr_script}")
    se      <- createSE(
        sample       = "${sample_id}",
        gene_counts  = readRDS("${gene_counts}"),
        exon_counts  = readRDS("${exon_counts}"),
        annotation   = readRDS("${intron_annotation}")
    )
    rpkm_se <- getGeneSE(se, rpkm.cutoff = 0.5)
    expressed <- getExpressedGenes(rpkm_se)
    saveRDS(expressed,    "${sample_id}_rpkm.rds")
    saveRDS(se,           "${sample_id}_all.rds")
    REOF
    """
}
