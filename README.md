# IPAseek

IPAseek is a bioinformatics pipeline to detect **intronic polyadenylation (IPA)** sites from RNA-seq datasets. It identifies locations within gene introns where premature polyadenylation occurs, quantifies IPA usage across samples, and produces a `SummarizedExperiment` object for downstream analysis.

> **Two execution modes are available:**
> - **Nextflow DSL2 pipeline** (recommended) — portable, reproducible, runs on any SLURM cluster with a single command. See [`nextflow/README.md`](nextflow/README.md).
> - **Original R/SLURM scripts** — the underlying R scripts described below, used directly on an HPC cluster.

---

## Quick Start

### Step 1 — One-time setup (run once per cluster)

```bash
bash setup.sh
```

This will ask for 4 things:

| Prompt | What to provide |
|--------|----------------|
| Project directory | Absolute path where IPAseek is installed |
| STAR index path | Path to your pre-built hg38 STAR genome index (~30 GB) |
| Annotation directory | Directory containing `hg38.cds.rds` for gene counting |
| SLURM account | Your cluster account name (leave blank if not required) |

### Step 2 — Run the pipeline

```bash
source ipaseek.env   # load configuration

nextflow run nextflow/main.nf \
    --data_table input_data_tables/data_table_test.txt \
    --atlas_name my_experiment \
    --outdir results
```

That's it. The pipeline reads your `data_table.txt`, runs all 3 stages, and produces a final `SummarizedExperiment` object in `results/`.

> **Alternative:** use the original R/SLURM scripts (see [Running the Original R/SLURM Pipeline](#running-the-original-rslurm-pipeline) below).

---

## Repository Structure

```
IPAseek/
├── nextflow/                      # ★ Nextflow DSL2 pipeline (recommended)
│   ├── main.nf                    # Entry point
│   ├── nextflow.config            # Pipeline configuration
│   ├── conf/                      # SLURM + resource configs
│   ├── modules/                   # 16 process definitions (stage1–3)
│   ├── workflows/                 # Stage sub-workflows
│   ├── assets/                    # Samplesheet template
│   └── README.md                  # Nextflow-specific documentation
├── 1_intron_preprocessing/        # Pre-computed hg38 intron annotation objects
│   ├── 1_flatten_genome/          # Flattened genome annotation (.rds)
│   ├── 2_annotation_object/       # CDS-level intron annotation objects (.rds)
│   └── 3_filtering_gobj/          # Filtered intron sets and filter files
├── 2_gene_preprocessing/          # Gene expression preprocessing
│   ├── 1_rnaseq_pipeline/         # STAR alignment + count scripts
│   ├── 2_bams/                    # BAM filtering scripts
│   ├── 3_gene_expression/         # Gene expression SE result objects
│   └── run_gene_exp.R             # Top-level orchestration script (Stage 2)
├── 3_ipa_run/
│   └── scripts/                   # IPA detection scripts (Stages 3 steps 1–11)
│       ├── 0_ipa_detect_run.r     # Top-level orchestration script (Stage 3)
│       ├── 1_ipa_detect.r         # Genomic coverage & intron retention
│       ├── 2_filtering.r          # Coverage-based filtering
│       ├── 3_pelt.r               # PELT changepoint detection
│       ├── 4_filter_cpts_de.R     # Changepoint filtering
│       ├── 5_merge.r              # Merge IPA calls across chromosomes
│       ├── 6_analyse_exon_structure.r  # TPM for new terminal exons
│       ├── 7_make_atlas.r         # Build IPA atlas per sample group
│       ├── 8_calculate_ipa_usage_combined.R  # IPA usage calculation
│       ├── 9_create_ipa_usage_se.R           # IPA usage SummarizedExperiment
│       └── master_workflow.R      # SLURM-aware wrapper (optional)
├── input_data_tables/             # Sample metadata tables for test datasets
│   ├── data_table_test.txt
│   └── data_table_test2.txt
├── environment.yml                # Conda environment (R + system tools)
├── CITATION.cff                   # Machine-readable citation
├── LICENSE                        # MIT License
├── pelt/                          # PELT output directory (created at runtime)
└── test_input/                    # Small test BAM/FASTQ files
```

---

## Requirements

### System

- [Nextflow](https://www.nextflow.io/) ≥ 23.04 and Java ≥ 11 *(for Nextflow mode)*
- A SLURM-based HPC cluster (job submission via `sbatch`)
- [STAR](https://github.com/alexdobin/STAR) aligner ≥ 2.7.10a
- [samtools](http://www.htslib.org/) ≥ 1.17
- A STAR genome index for the reference genome (e.g., hg38)

### R packages

Install all dependencies at once using conda:

```bash
conda env create -f environment.yml
conda activate ipaseek
```

Or install manually in R:

```r
install.packages(c("data.table", "gtools", "dplyr", "tidyr", "tidyverse", "ggplot2", "scales", "changepoint", "slurmR"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GenomicAlignments",
  "GenomicFeatures",
  "GenomicRanges",
  "SummarizedExperiment",
  "DESeq2"
))
```

---

## Running the Original R/SLURM Pipeline

> If you prefer the Nextflow pipeline, see [`nextflow/README.md`](nextflow/README.md).

Set the project directory and data table variables at the top of each orchestration script before running.

### Input Data Table Format

All pipeline stages are driven by a tab-delimited sample metadata table:

| Column | Description |
|--------|-------------|
| `FILE_PATH` | Directory where BAM files are (or will be) stored |
| `UNIQUE_ID` | SRR accession or other unique identifier for the sample |
| `NAME` | Sample name used throughout the pipeline |
| `GENOME` | Reference genome name (e.g., `hg38`) |
| `FASTQ_FILE` | Path(s) to FASTQ file(s); paired-end files separated by `::` |
| `CONDITION` | Experimental condition label |

Example (`input_data_tables/data_table_test.txt`):

```
FILE_PATH	UNIQUE_ID	NAME	GENOME	FASTQ_FILE	CONDITION
/path/to/bams/uniq_bams	SRR6830273	CLL7	hg38	/path/to/SRR6830273_1.fastq.gz::/path/to/SRR6830273_2.fastq.gz	YES
/path/to/bams/uniq_bams	SRR6830274	CLL11	hg38	/path/to/SRR6830274_1.fastq.gz::/path/to/SRR6830274_2.fastq.gz	YES
```

### Stage 1: Intron Preprocessing (pre-computed)

Pre-computed hg38 intron annotation objects are already provided under `1_intron_preprocessing/`. No action is required for hg38 datasets. For other genomes, equivalent annotation `.rds` files must be prepared and placed in the corresponding subdirectories.

### Stage 2: Gene Expression Preprocessing (`run_gene_exp.R`)

```r
source("2_gene_preprocessing/run_gene_exp.R")
```

| Step | Function | Description |
|------|----------|-------------|
| 1 | `pushSTAR()` | Aligns FASTQ files to the reference genome using STAR |
| 2 | `filt_bam()` | Filters BAM files to retain uniquely mapped reads |
| 3 | `runCDSCounts()` | Counts reads over CDS regions for each sample |
| 4 | `createSE()` / `getGeneSE()` / `getExpressedGenes()` | Builds `SummarizedExperiment` objects for raw counts and RPKM; saves per-sample `.rds` files |
| 5 | `retrieve_geneexpr_data()` / `geneexpr_se()` | Retrieves and consolidates expression values from the SE objects |

**Key output:** Per-sample gene expression `.rds` files saved to
`2_gene_preprocessing/3_gene_expression/results/<sample_name>/gene_expression/`.

### Stage 3: IPA Detection (`3_ipa_run/scripts/0_ipa_detect_run.r`)

```r
source("3_ipa_run/scripts/0_ipa_detect_run.r")
```

| Step | Function | Description |
|------|----------|-------------|
| 1 | `ipa_detect()` | Calculates per-base genomic coverage over introns; submits per-sample SLURM jobs |
| 2 | `retrieve_intronreten_data()` | Collects intron retention data from completed jobs |
| 3 | `intronret_se()` | Builds a `SummarizedExperiment` for intron retention across samples |
| 4 | `filtering()` | Filters introns requiring ≥5 reads over a continuous ≥100 bp stretch; submits per-chromosome SLURM jobs |
| 5 | `pelt()` | Runs PELT changepoint detection on coverage signals; submits per-chromosome SLURM jobs |
| 6 | `filter_changepoints()` | Filters changepoints to identify high-confidence IPA events; submits per-sample SLURM jobs |
| 7 | `merge_cpts()` | Merges per-chromosome IPA calls into per-sample genome-wide results |
| 8 | `filter_te_run()` | Analyses terminal exon structure and computes TPM for new terminal exons |
| 9 | `make_atlas_run()` | Constructs an IPA site atlas for each sample group |
| 10 | `calc_ipa_usage()` | Quantifies IPA usage (reads, RPKM) for each site in every sample |
| 11 | `ipa_usage_se()` | Assembles all assays into a final `SummarizedExperiment` object |

**Key outputs** (written under `pelt/results/<atlas_name>/`):

- `<atlas_name>_full_ipa_atlas_conf.RDS` — confident IPA site atlas (GRanges)
- `<atlas_name>_<sample>_ipa_usage_atlas.csv` — per-sample IPA usage table
- `<atlas_name>_ipa_usage_se.RDS` — multi-assay `SummarizedExperiment` with IPA usage, read counts, and RPKM across all samples

---

## Test Run with Example Data

The `input_data_tables/` directory contains two small test metadata tables:

- `data_table_test.txt` — used for alignment and BAM generation (Stage 2)
- `data_table_test2.txt` — alternative test dataset

For Nextflow, edit `nextflow/assets/samplesheet_template.csv` with your paths and run:

```bash
nextflow run nextflow/main.nf -profile slurm \
    --samplesheet nextflow/assets/samplesheet_template.csv \
    --star_index /path/to/STAR_index \
    --intron_annotation 1_intron_preprocessing/3_filtering_gobj/rnhg38_filtered_introns_cds.rds \
    --outdir results
```

---

## Notes

- IPAseek uses SLURM (`sbatch`) to parallelize compute-intensive steps across chromosomes and samples.
- The `master_workflow.R` script provides an alternative orchestration layer that polls SLURM job status between steps using `sacct`.
- All pre-computed annotation objects shipped in `1_intron_preprocessing/` are for **hg38**. Analyses on other genome builds require regenerating these objects.
- The `pelt/` directory is created automatically at runtime and holds all intermediate and final results.

---

## Citation

If you use IPAseek in your research, please cite:

> Rashmi *et al.* (2026). *Dynamics of intronic polyadenylation in the hematopoietic lineage and its regulation by DNA methylation*. **Genome Research**, 36(6):1250. https://genome.cshlp.org/content/36/6/1250

A machine-readable citation is available via the **"Cite this repository"** button on GitHub (powered by [`CITATION.cff`](CITATION.cff)).
