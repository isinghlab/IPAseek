# IPAseek — Nextflow DSL2 Pipeline

This directory contains a **Nextflow DSL2** version of the IPAseek intronic-polyadenylation (IPA) detection pipeline. It faithfully reproduces all three stages and sixteen processes of the original SLURM-based R pipeline.

---

## Prerequisites

| Tool | Minimum version |
|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.04 |
| Java | ≥ 11 |
| SLURM | any recent version |
| [STAR](https://github.com/alexdobin/STAR) | ≥ 2.7.10a |
| [samtools](http://www.htslib.org/) | ≥ 1.17 |
| R | ≥ 4.3.0 |
| R packages | see `../environment.yml` |
| Conda / Mamba | optional but recommended |

Install the full R + tool environment with:

```bash
conda env create -f ../environment.yml
conda activate ipaseek
```

---

## Quick start

```bash
nextflow run nextflow/main.nf \
    -profile slurm \
    --samplesheet assets/samplesheet_template.csv \
    --star_index /path/to/STAR_genome_index \
    --intron_annotation 1_intron_preprocessing/3_filtering_gobj/rnhg38_filtered_introns_cds.rds \
    --outdir results
```

### Re-running Stage 1 (intron filtering from scratch)

Omit `--intron_annotation` and supply the annotation directory instead:

```bash
nextflow run nextflow/main.nf \
    -profile slurm \
    --samplesheet assets/samplesheet_template.csv \
    --star_index /path/to/STAR_genome_index \
    --annotation_dir 1_intron_preprocessing/3_filtering_gobj \
    --outdir results
```

---

## Samplesheet format

Create a CSV file following the template in `assets/samplesheet_template.csv`:

```
sample,condition,fastq_1,fastq_2,genome
CLL7,treatment,/path/to/SRR6830273_1.fastq.gz,/path/to/SRR6830273_2.fastq.gz,hg38
CLL11,control,/path/to/SRR6830274_1.fastq.gz,/path/to/SRR6830274_2.fastq.gz,hg38
```

| Column | Description |
|---|---|
| `sample` | Unique sample identifier |
| `condition` | Experimental condition (e.g. treatment / control) |
| `fastq_1` | Absolute path to Read 1 FASTQ (can be gzipped) |
| `fastq_2` | Absolute path to Read 2 FASTQ (can be gzipped) |
| `genome` | Genome build (currently only `hg38` is supported) |

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | *required* | Path to input CSV samplesheet |
| `--genome` | `hg38` | Genome build name |
| `--star_index` | *required* | Path to STAR genome index directory |
| `--intron_annotation` | `null` | Path to pre-computed `rnhg38_filtered_introns_cds.rds`; if omitted, Stage 1 re-runs |
| `--annotation_dir` | `null` | Path to annotation directory (used when `--intron_annotation` is not provided) |
| `--atlas_name` | `ipaseek_atlas` | Prefix for IPA atlas output files |
| `--outdir` | `results` | Directory for all pipeline outputs |
| `--chromosomes` | chr1–22, chrX, chrY | List of chromosomes to process in Stage 3 |

---

## Pipeline stages

### Stage 1 — Intron preprocessing

| Process | Script |
|---|---|
| `FILTER_INTRONS` | `1_intron_preprocessing/3_filtering_gobj/scripts/1_filtering_gobj.r` |

Produces `rnhg38_filtered_introns_cds.rds` — the core intron annotation used by all downstream stages.

### Stage 2 — Gene expression preprocessing

| Process | Description |
|---|---|
| `STAR_ALIGN` | STAR 2-pass alignment (paired-end, sorted BAM) |
| `FILTER_BAM` | Retain uniquely mapping reads (NH:i:1) |
| `COUNT_GENES` | Count reads over genes and exons |
| `GENE_EXPRESSION_SE` | Build per-sample SummarizedExperiment; filter by RPKM ≥ 0.5 |
| `MERGE_GENE_EXPR` | Merge all samples into cohort-level SE objects |

### Stage 3 — IPA detection

| Process | Description |
|---|---|
| `IPA_DETECT` | Compute per-chr coverage and intron retention |
| `COLLECT_RETENTION` | Aggregate retention files across chromosomes per sample |
| `INTRON_RETENTION_SE` | Build intron-retention SummarizedExperiment |
| `FILTER_COVERAGE` | Filter coverage by retention threshold |
| `PELT_CHANGEPOINT` | PELT changepoint detection (pen=c(100,10000), minseglen=200) |
| `FILTER_CHANGEPOINTS` | Filter changepoints by expression level |
| `MERGE_CHANGEPOINTS` | Merge per-chromosome changepoint results per sample |
| `FILTER_TERMINAL_EXONS` | Analyse exon structure; filter terminal exons |
| `MAKE_ATLAS` | Build cohort-level IPA atlas |
| `CALC_IPA_USAGE` | Calculate per-sample IPA usage from atlas |
| `IPA_USAGE_SE` | Build final IPA usage SummarizedExperiment |

---

## Output structure

```
results/
├── pipeline_info/          # Nextflow execution timeline, report, trace, DAG
├── stage1/
│   └── filter_introns/     # rnhg38_filtered_introns_cds.rds
├── stage2/
│   ├── star_align/         # Per-sample sorted BAM + BAI
│   ├── filter_bam/         # Per-sample uniquely-mapped BAM + BAI
│   ├── count_genes/        # Per-sample gene/exon count RDS
│   ├── gene_expression_se/ # Per-sample RPKM and full SE RDS
│   └── merge_gene_expr/    # se_gene_expr_count.RDS, se_gene_expr_rpkm.RDS
└── stage3/
    ├── ipa_detect/         # Per-sample × chr coverage + retention RDS
    ├── collect_retention/  # Per-sample aggregated retention
    ├── intron_retention_se/ # se_intron_retention.RDS
    ├── filter_coverage/    # Per-sample × chr filtered coverage RDS
    ├── pelt_changepoint/   # Per-sample × chr PELT results RDS
    ├── filter_changepoints/ # Per-sample filtered changepoint CSVs
    ├── merge_changepoints/ # Per-sample merged changepoint CSVs
    ├── filter_terminal_exons/ # Per-sample terminal exon CSVs
    ├── make_atlas/         # Cohort IPA atlas RDS + CSV
    ├── calc_ipa_usage/     # Per-sample IPA usage CSV
    └── ipa_usage_se/       # Final IPA usage SE RDS
```

---

## SLURM configuration

Edit `conf/slurm.config` to set your cluster account and queue:

```groovy
process {
    executor       = 'slurm'
    clusterOptions = '--account=YOUR_ACCOUNT'
    queue          = 'normal'
}
```

Resource allocations (CPUs, memory, walltime) are defined in `conf/resources.config` and mirror the original SLURM submission scripts.

---

## Tests

See `.nf-test/README.md` for instructions on adding nf-test unit and integration tests.

---

## Citation

If you use IPAseek, please cite:

> Rashmi *et al.* (2026). *Cancer-associated dynamics and potential regulators of intronic polyadenylation revealed by IPAFinder using standard RNA-seq data*. **Genome Research**, 36(6):1250. https://doi.org/10.1101/gr.278518.123
