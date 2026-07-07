process MERGE_GENE_EXPR {
    tag 'merge_gene_expr'
    publishDir "${params.outdir}/stage2/merge_gene_expr", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    path intron_annotation
    path rpkm_files

    output:
    path "se_gene_expr_count.RDS", emit: se_count
    path "se_gene_expr_rpkm.RDS",  emit: se

    script:
    def expr_script = "${projectDir}/../../2_gene_preprocessing/1_rnaseq_pipeline/expressed_genes_data_retrieved.R"
    """
    Rscript - <<'REOF'
    source("${expr_script}")
    rpkm_files   <- list.files(".", pattern = "_rpkm\\\\.rds\$", full.names = TRUE)
    annotation   <- readRDS("${intron_annotation}")
    gene_data    <- retrieve_geneexpr_data(rpkm_files, annotation)
    se_count     <- geneexpr_se(gene_data, type = "count")
    se_rpkm      <- geneexpr_se(gene_data, type = "rpkm")
    saveRDS(se_count, "se_gene_expr_count.RDS")
    saveRDS(se_rpkm,  "se_gene_expr_rpkm.RDS")
    REOF
    """
}
