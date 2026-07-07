# IPAseek

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2Fgr.278518.123-blue)](https://genome.cshlp.org/content/36/6/1250)
[![R ≥ 4.3](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3.svg?logo=r)](https://www.r-project.org/)

IPAseek is a bioinformatics pipeline to detect **intronic polyadenylation (IPA)** sites from RNA-seq datasets. It identifies locations within gene introns where premature polyadenylation occurs, quantifies IPA usage across samples, and produces a `SummarizedExperiment` object for downstream analysis.

---

## Quick Start

### 1 — Install dependencies (once)

```bash
conda env create -f environment.yml
conda activate ipaseek
```

### 2 — Configure for your cluster (once)

```bash
bash setup.sh
```

`setup.sh` asks for **4 things** and writes an `ipaseek.env` config file:

| Prompt | What to provide |
|--------|----------------|
| **Project directory** | Absolute path where IPAseek is cloned |
| **STAR index path** | Path to your pre-built hg38 STAR genome index (~30 GB) |
| **Annotation directory** | Directory containing `hg38.cds.rds` for gene counting |
| **SLURM account** | Your cluster account name (leave blank if not required) |

> 💡 `setup.sh` auto-detects your SLURM account and checks all required tools are on `PATH`.

### 3 — Run the pipeline

```bash
source ipaseek.env

nextflow run nextflow/main.nf \
    --data_table input_data_tables/data_table_test.txt \
    --atlas_name my_experiment \
    --outdir results
```

**That's it.** Point to your `data_table.txt`, the pipeline runs all 3 stages end-to-end and produces a final `SummarizedExperiment` object in `results/`. The intermediate `_uniq` data table is generated automatically — no manual steps.

---

## Input: `data_table.txt`

The only file you need to prepare. A tab-delimited table with one row per sample:

| Column | Description | Example |
|--------|-------------|---------|
| `FILE_PATH` | Directory where BAM files will be stored | `/path/to/bams` |
| `UNIQUE_ID` | SRR accession or unique sample ID | `SRR6830273` |
| `NAME` | Sample name used throughout the pipeline | `CLL7` |
| `GENOME` | Reference genome (`hg38` supported) | `hg38` |
| `FASTQ_FILE` | FASTQ path(s); paired-end separated by `::` | `/path/R1.fastq.gz::/path/R2.fastq.gz` |
| `CONDITION` | Experimental condition label | `YES` / `NO` |

A template is at [`input_data_tables/data_table_test.txt`](input_data_tables/data_table_test.txt).

---

## What the pipeline does

```
data_table.txt
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│ Stage 1 — Intron Preprocessing                                  │
│   Pre-computed hg38 annotation (ships with the repo — no action)│
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│ Stage 2 — Gene Expression Preprocessing            per sample   │
│   STAR align → filter unique BAMs → count genes →              │
│   SummarizedExperiment → merge across samples                   │
│   (uniq data table auto-generated after BAM filtering)          │
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│ Stage 3 — IPA Detection                        per sample × chr │
│   Coverage → intron retention → PELT changepoints →             │
│   filter → merge → terminal exon analysis →                     │
│   IPA atlas → usage quantification → final SE                   │
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
results/<atlas_name>_ipa_usage_se.RDS   ← final output
```

### Stage 2 steps

| Step | Function | Description |
|------|----------|-------------|
| 1 | `pushSTAR()` | Aligns FASTQ files using STAR (paired-end, sorted BAM) |
| 2 | `filt_bam()` | Filters BAMs to uniquely mapping reads (NH:i:1) |
| — | *(auto)* | Generates `data_table_uniq.txt` with updated paths + library sizes |
| 3 | `runCDSCounts()` | Counts reads over CDS regions per sample |
| 4 | `createSE()` / `getExpressedGenes()` | Builds per-sample SummarizedExperiment (RPKM ≥ 0.5 filter) |
| 5 | `geneexpr_se()` | Merges all samples into cohort-level SE objects |

### Stage 3 steps

| Step | Function | Description |
|------|----------|-------------|
| 1 | `ipa_detect()` | Per-base coverage over introns; per-sample × chr SLURM jobs |
| 2 | `retrieve_intronreten_data()` | Collects retention data across chromosomes |
| 3 | `intronret_se()` | Intron retention SummarizedExperiment |
| 4 | `filtering()` | Coverage filter: ≥5 reads over ≥100 bp continuous stretch |
| 5 | `pelt()` | PELT changepoint detection (pen=c(100,10000), minseglen=200) |
| 6 | `filter_changepoints()` | High-confidence IPA event selection |
| 7 | `merge_cpts()` | Per-sample genome-wide merge |
| 8 | `filter_te_run()` | Terminal exon structure + TPM |
| 9 | `make_atlas_run()` | Cohort-level IPA site atlas |
| 10 | `calc_ipa_usage()` | IPA usage quantification (reads, RPKM) per site per sample |
| 11 | `ipa_usage_se()` | Final multi-assay SummarizedExperiment |

---

## Outputs

All outputs are written to `--outdir` (default: `results/`):

```
results/
├── stage2/
│   ├── star_align/           # Per-sample sorted BAM + BAI
│   ├── filter_bam/           # Uniquely-mapped BAM + BAI
│   ├── count_genes/          # Per-sample gene/exon count RDS
│   ├── gene_expression_se/   # Per-sample RPKM RDS
│   └── merge_gene_expr/      # se_gene_expr_count.RDS, se_gene_expr_rpkm.RDS
└── stage3/
    ├── ...intermediate files...
    ├── make_atlas/           # <atlas_name>_full_ipa_atlas_conf.RDS + .csv
    ├── calc_ipa_usage/       # Per-sample IPA usage CSV
    └── ipa_usage_se/         # ★ <atlas_name>_ipa_usage_se.RDS  ← final output
```

**Key final outputs:**
- `<atlas_name>_full_ipa_atlas_conf.RDS` — confident IPA site atlas (GRanges)
- `<atlas_name>_<sample>_ipa_usage_atlas.csv` — per-sample IPA usage table
- `<atlas_name>_ipa_usage_se.RDS` — multi-assay SummarizedExperiment (IPA usage, read counts, RPKM)

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.04 | + Java ≥ 11 |
| [STAR](https://github.com/alexdobin/STAR) | ≥ 2.7.10a | For RNA-seq alignment |
| [samtools](http://www.htslib.org/) | ≥ 1.17 | For BAM processing |
| R | ≥ 4.3.0 | Via conda (recommended) |
| SLURM | any | HPC job scheduler |

Install all R packages and tools at once:
```bash
conda env create -f environment.yml
conda activate ipaseek
```

---

## Repository Structure

```
IPAseek/
├── setup.sh                       # ★ Start here — one-time setup
├── ipaseek.env                    # Generated by setup.sh (not committed)
├── nextflow/                      # Nextflow DSL2 pipeline
│   ├── main.nf                    # Entry point (reads data_table.txt natively)
│   ├── nextflow.config            # Parameters + env var integration
│   ├── conf/
│   │   ├── slurm.config           # SLURM executor (account auto-read from env)
│   │   ├── resources.config       # Per-process CPU/memory/time
│   │   └── test.config            # Smoke test profile (chr22 only)
│   ├── modules/                   # 16 process definitions
│   ├── workflows/                 # Stage sub-workflows
│   ├── assets/samplesheet_template.csv
│   └── README.md                  # Nextflow-specific docs
├── 1_intron_preprocessing/        # Pre-computed hg38 annotation (ships with repo)
├── 2_gene_preprocessing/          # R scripts: STAR, BAM filtering, gene SE
│   └── run_gene_exp.R             # Stage 2 orchestration (accepts CLI args)
├── 3_ipa_run/scripts/             # R scripts: IPA detection (11 steps)
│   └── 0_ipa_detect_run.r         # Stage 3 orchestration (accepts CLI args)
├── input_data_tables/             # data_table_test.txt template
├── test_input/                    # Small test BAM/FASTQ files
├── environment.yml                # Conda environment
├── CITATION.cff                   # Machine-readable citation
└── LICENSE                        # MIT
```

---

## Running with Original R/SLURM Scripts (Advanced)

If you prefer to run the R scripts directly instead of Nextflow:

```bash
source ipaseek.env

# Stage 2: Gene expression
Rscript 2_gene_preprocessing/run_gene_exp.R \
    --project_dir $IPASEEK_PROJECT_DIR \
    --data_table  input_data_tables/data_table_test.txt \
    --atlas_name  my_experiment

# Stage 3: IPA detection (uses auto-generated _uniq table from Stage 2)
Rscript 3_ipa_run/scripts/0_ipa_detect_run.r \
    --project_dir $IPASEEK_PROJECT_DIR \
    --data_table  input_data_tables/data_table_test_uniq.txt \
    --atlas_name  my_experiment
```

---

## Smoke Test

Verify your setup with the bundled test data (chr22 only, completes in minutes):

```bash
source ipaseek.env

nextflow run nextflow/main.nf \
    -profile test,slurm \
    --outdir test_results
```

---

## Notes

- All pre-computed annotation objects in `1_intron_preprocessing/` are for **hg38**. Other genome builds require regenerating these objects.
- The `pelt/` directory (R/SLURM mode) and `results/` (Nextflow mode) are created automatically at runtime.
- SLURM parallelises compute-intensive steps across chromosomes and samples automatically.

---

## Citation

If you use IPAseek in your research, please cite:

> Rashmi *et al.* (2026). *Cancer-associated dynamics and potential regulators of intronic polyadenylation revealed by IPAFinder using standard RNA-seq data*. **Genome Research**, 36(6):1250. https://genome.cshlp.org/content/36/6/1250

A machine-readable citation is available via the **"Cite this repository"** button on GitHub (powered by [`CITATION.cff`](CITATION.cff)).
