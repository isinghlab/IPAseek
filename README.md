# IPAseek

IPAseek is a bioinformatics pipeline to detect **intronic polyadenylation (IPA)** sites from RNA-seq datasets.

This repository contains:

- The original IPAseek scripts under `IPAseek/`
- A SLURM-based workflow layer under `workflow/` to make the pipeline easy to run on HPC clusters
- Example input tables under `input_data_tables/`
- Example configuration files under `examples/example_config/`

## Requirements

- A SLURM-based HPC cluster (for job submission via `sbatch`)
- R (version X.Y.Z+) and required packages
- Optional: Conda for environment management

## Quickstart

1. **Clone the repository**

```bash
git clone https://github.com/richarashmivats/IPAseek.git
cd IPAseek
```

2. **(Optional) Create a conda environment**

```bash
conda env create -f environment/environment.yml
conda activate ipaseek
```

3. **Copy and edit example config files**

```bash
cp examples/example_config/config.yaml config/config.yaml
cp examples/example_config/samples.csv config/samples.csv
```

Edit:

- `config/config.yaml` to point to your genome FASTA, GTF, and to set SLURM parameters (account, partition, etc.).
- `config/samples.csv` to list your RNA-seq samples and BAM file paths.

4. **Submit the pipeline**

```bash
bash workflow/run_ipaseek_slurm.sh -c config/config.yaml -s config/samples.csv
```

This will submit three dependent SLURM jobs:

1. Intron preprocessing
2. Gene expression preprocessing
3. IPA detection

You can monitor them with `squeue` as usual on your cluster.

## Test run with example data

The `input_data_tables/` directory includes small test tables such as:

- `data_table_test.txt`
- `data_table_test2.txt`

The example `config.yaml` uses `data_table_test_uniq.txt` so you can run a quick test after editing only genome and SLURM settings.

```bash
bash workflow/run_ipaseek_slurm.sh \
  -c examples/example_config/config.yaml \
  -s examples/example_config/samples.csv
```

## Notes

- IPAseek currently supports SLURM clusters and submits jobs with `sbatch`. Running on non-SLURM systems will require manually adapting the workflow.
- The high-level workflow code is under `workflow/`. The core IPAseek logic remains in the original `IPAseek/` subdirectories.
