include { IPA_DETECT           } from '../modules/stage3/ipa_detect'
include { COLLECT_RETENTION    } from '../modules/stage3/collect_retention'
include { INTRON_RETENTION_SE  } from '../modules/stage3/intron_retention_se'
include { FILTER_COVERAGE      } from '../modules/stage3/filter_coverage'
include { PELT_CHANGEPOINT     } from '../modules/stage3/pelt_changepoint'
include { FILTER_CHANGEPOINTS  } from '../modules/stage3/filter_changepoints'
include { MERGE_CHANGEPOINTS   } from '../modules/stage3/merge_changepoints'
include { FILTER_TERMINAL_EXONS } from '../modules/stage3/filter_terminal_exons'
include { MAKE_ATLAS           } from '../modules/stage3/make_atlas'
include { CALC_IPA_USAGE       } from '../modules/stage3/calc_ipa_usage'
include { IPA_USAGE_SE         } from '../modules/stage3/ipa_usage_se'

workflow STAGE3 {
    take:
    uniq_bams_ch      // channel: [sample_id, condition, bam, bai]
    rpkm_ch           // channel: [sample_id, condition, rpkm_rds]
    intron_annotation // path to intron annotation RDS
    atlas_name        // string: atlas name

    main:
    // Create a per-sample × per-chromosome channel
    chromosomes_ch = Channel.from(params.chromosomes)
    bam_chr_ch = uniq_bams_ch.combine(chromosomes_ch)

    // Flatten rpkm_ch to just [sample_id, rpkm_rds] for joining
    rpkm_flat_ch = rpkm_ch.map { sample_id, condition, rpkm -> tuple(sample_id, rpkm) }

    // Stage 3a: IPA coverage and retention per sample × chromosome
    IPA_DETECT(bam_chr_ch, intron_annotation, rpkm_flat_ch.map { it[1] })

    // Collect per-chromosome retention files per sample
    COLLECT_RETENTION(
        IPA_DETECT.out.retention
            .map { sample_id, chr, f -> tuple(sample_id, f) }
            .groupTuple(by: 0)
    )

    // Build intron-retention SummarizedExperiment across all samples
    INTRON_RETENTION_SE(
        COLLECT_RETENTION.out.retention_data
            .map { sample_id, f -> f }
            .collect()
    )

    // Filter coverage per sample × chromosome (join retention summary)
    filter_cov_input_ch = IPA_DETECT.out.coverage
        .join(IPA_DETECT.out.retention, by: [0, 1])
        .map { sample_id, chr, cov_rds, retain_rds ->
            tuple(sample_id, chr, cov_rds, retain_rds)
        }
        .combine(
            COLLECT_RETENTION.out.retention_all.map { sample_id, f -> tuple(sample_id, f) },
            by: 0
        )
        .map { sample_id, chr, cov_rds, retain_rds, retention_summary ->
            tuple(sample_id, chr, cov_rds, retain_rds, retention_summary)
        }

    FILTER_COVERAGE(filter_cov_input_ch)

    // PELT changepoint detection per sample × chromosome
    pelt_input_ch = FILTER_COVERAGE.out.filtered
        .join(IPA_DETECT.out.coverage, by: [0, 1])

    PELT_CHANGEPOINT(pelt_input_ch)

    // Filter changepoints per sample (join grouped PELT results with BAM info)
    filter_cpts_input_ch = PELT_CHANGEPOINT.out.pelt
        .map { sample_id, chr, f -> tuple(sample_id, f) }
        .groupTuple(by: 0)
        .join(uniq_bams_ch.map { sample_id, condition, bam, bai -> tuple(sample_id, condition, bam, bai) }, by: 0)
        .map { sample_id, pelt_files, condition, bam, bai ->
            tuple(sample_id, condition, pelt_files, bam, bai)
        }

    rpkm_for_filter_ch = rpkm_flat_ch.map { sample_id, rpkm -> rpkm }

    FILTER_CHANGEPOINTS(filter_cpts_input_ch, intron_annotation, rpkm_for_filter_ch)

    // Merge changepoints per sample
    MERGE_CHANGEPOINTS(FILTER_CHANGEPOINTS.out.cpts)

    // Filter terminal exons per sample (join with BAM info)
    te_input_ch = MERGE_CHANGEPOINTS.out.merged
        .join(uniq_bams_ch.map { sample_id, condition, bam, bai -> tuple(sample_id, bam, bai) }, by: 0)

    FILTER_TERMINAL_EXONS(te_input_ch, intron_annotation)

    // Build IPA atlas from all samples
    MAKE_ATLAS(
        FILTER_TERMINAL_EXONS.out.te.map { sample_id, f -> f }.collect(),
        intron_annotation,
        atlas_name
    )

    // Calculate IPA usage per sample
    CALC_IPA_USAGE(uniq_bams_ch, MAKE_ATLAS.out.atlas_rds, intron_annotation)

    // Build final IPA usage SummarizedExperiment
    IPA_USAGE_SE(
        CALC_IPA_USAGE.out.usage.map { sample_id, f -> f }.collect(),
        atlas_name
    )

    emit:
    ipa_se = IPA_USAGE_SE.out.se
}
