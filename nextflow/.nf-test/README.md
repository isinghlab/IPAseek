# nf-test placeholder

This directory is reserved for [nf-test](https://www.nf-test.com/) unit and integration tests for the IPAseek Nextflow pipeline.

## How to add tests

1. **Install nf-test** (requires Java 11+ and Nextflow ≥23.04):
   ```bash
   curl -fsSL https://get.nf-test.com | bash
   ```

2. **Initialise nf-test in this directory** (run from `nextflow/`):
   ```bash
   nf-test init
   ```

3. **Generate a test for a module**, e.g. `FILTER_BAM`:
   ```bash
   nf-test generate process modules/stage2/filter_bam.nf
   ```

4. **Run all tests**:
   ```bash
   nf-test test
   ```

## Suggested tests to write

| Process / Workflow | Test description |
|---|---|
| `FILTER_INTRONS` | Check that `rnhg38_filtered_introns_cds.rds` is produced |
| `STAR_ALIGN` | Check BAM + BAI outputs exist for a tiny FASTQ pair |
| `FILTER_BAM` | Verify only NH:i:1 reads remain in output BAM |
| `COUNT_GENES` | Check `*_gene_counts.rds` and `*_exon_counts.rds` outputs |
| `IPA_DETECT` | Check per-chromosome coverage / retention RDS files |
| `PELT_CHANGEPOINT` | Check `ipa_final_*_filtered_*.RDS` output |
| `MAKE_ATLAS` | Verify atlas RDS and CSV are produced |
| Full `STAGE2` workflow | End-to-end with two test samples |
| Full `STAGE3` workflow | End-to-end with pre-computed atlas inputs |
