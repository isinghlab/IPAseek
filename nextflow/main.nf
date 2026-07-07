#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { STAGE1 } from './workflows/stage1'
include { STAGE2 } from './workflows/stage2'
include { STAGE3 } from './workflows/stage3'

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
    STAGE3(
        STAGE2.out.uniq_bams,
        STAGE2.out.rpkm_per_sample,
        intron_annotation,
        params.atlas_name
    )
}
