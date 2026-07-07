include { STAR_ALIGN        } from '../modules/stage2/star_align'
include { FILTER_BAM        } from '../modules/stage2/filter_bam'
include { COUNT_GENES       } from '../modules/stage2/count_genes'
include { GENE_EXPRESSION_SE } from '../modules/stage2/gene_expression_se'
include { MERGE_GENE_EXPR   } from '../modules/stage2/merge_gene_expr'

workflow STAGE2 {
    take:
    reads_ch          // channel: [sample_id, condition, fastq_1, fastq_2]
    star_index        // path to STAR genome index
    intron_annotation // path to intron annotation RDS

    main:
    STAR_ALIGN(reads_ch, star_index)
    FILTER_BAM(STAR_ALIGN.out.bam)
    COUNT_GENES(FILTER_BAM.out.uniq_bam, intron_annotation)
    GENE_EXPRESSION_SE(COUNT_GENES.out.counts, intron_annotation)
    MERGE_GENE_EXPR(intron_annotation, GENE_EXPRESSION_SE.out.rpkm.map { it[2] }.collect())

    emit:
    uniq_bams        = FILTER_BAM.out.uniq_bam
    rpkm_per_sample  = GENE_EXPRESSION_SE.out.rpkm
    gene_expr_se     = MERGE_GENE_EXPR.out.se
}
