#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { STAGE1 } from './workflows/stage1'
include { STAGE2 } from './workflows/stage2'
include { STAGE3 } from './workflows/stage3'

<<<<<<< HEAD
// Validate required params
def checkParams() {
    if (!params.data_table && !params.samplesheet) {
        error "Please provide either --data_table (IPAseek TSV) or --samplesheet (CSV)"
    }
    if (!params.star_index) {
        error "Please provide --star_index or run setup.sh to set IPASEEK_STAR_INDEX"
    }
}

workflow {
    checkParams()

    // ── Parse input: native data_table.txt or CSV samplesheet ──────────────
    if (params.data_table) {
        // Native IPAseek TSV format: FILE_PATH, UNIQUE_ID, NAME, GENOME, FASTQ_FILE, CONDITION
        reads_ch = Channel
            .fromPath(params.data_table)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def fastqs = row.FASTQ_FILE.tokenize('::')
                tuple(
                    row.NAME,
                    row.CONDITION,
                    file(fastqs[0]),
                    fastqs.size() > 1 ? file(fastqs[1]) : file(fastqs[0])
                )
            }
    } else {
        // Standard CSV samplesheet
        reads_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, row.condition, file(row.fastq_1), file(row.fastq_2)) }
    }

    // ── Stage 1: Intron annotation ──────────────────────────────────────────
    if (params.intron_annotation && file(params.intron_annotation).exists()) {
        intron_annotation = Channel.fromPath(params.intron_annotation)
    } else {
        STAGE1(Channel.fromPath("${projectDir}/../1_intron_preprocessing/3_filtering_gobj"))
        intron_annotation = STAGE1.out.intron_annotation
    }

    // ── Stage 2: Gene expression preprocessing ─────────────────────────────
    STAGE2(reads_ch, params.star_index, intron_annotation)

    // ── Stage 3: IPA detection ─────────────────────────────────────────────
=======
workflow {

    // ---------------------------------------------------------------------------
    // Parse samplesheet
    // ---------------------------------------------------------------------------
    reads_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.condition, file(row.fastq_1), file(row.fastq_2)) }

    // ---------------------------------------------------------------------------
    // Stage 1: Intron preprocessing
    // Use pre-computed annotation (params.intron_annotation) when provided,
    // otherwise re-run filtering via STAGE1.
    // ---------------------------------------------------------------------------
    if (params.intron_annotation) {
        intron_annotation = Channel.fromPath(params.intron_annotation)
    } else {
        STAGE1(Channel.fromPath(params.annotation_dir))
        intron_annotation = STAGE1.out.intron_annotation
    }

    // ---------------------------------------------------------------------------
    // Stage 2: Alignment and gene expression
    // ---------------------------------------------------------------------------
    STAGE2(reads_ch, params.star_index, intron_annotation)

    // ---------------------------------------------------------------------------
    // Stage 3: IPA detection
    // ---------------------------------------------------------------------------
>>>>>>> origin/main
    STAGE3(
        STAGE2.out.uniq_bams,
        STAGE2.out.rpkm_per_sample,
        intron_annotation,
        params.atlas_name
    )
}
