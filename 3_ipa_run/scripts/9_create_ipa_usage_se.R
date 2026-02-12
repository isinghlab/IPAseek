
library(GenomicAlignments)  
library(GenomicFeatures)    
library(dplyr)              
library(tidyr)              
library(data.table)         
library(gtools)             

ipa_usage_se <- function(input.data.path, wd, atlas_name) {

  # input.data.path <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/input_data_tables/data_table_test2.txt"
  # wd <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"
  # atlas_name <- "test2"
  # 
  ## Read the input sample metadata table (tab-delimited)
  data.input <- read.delim(input.data.path, sep="\t", header=TRUE)

  ## Load the confident IPA atlas (GRanges object) for this atlas
  ipa_atlas <- readRDS(paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_full_ipa_atlas_conf.RDS"))

  # Split and flatten the atlas to get individual IPA events
  ipa_atlas <- split(ipa_atlas, as.factor(ipa_atlas))
  ipa_atlas <- GRangesList(ipa_atlas)
  ipa_atlas <- unlist(ipa_atlas)

  ## Get the sample names included in this atlas
  sampleNames <- as.character(data.input$NAME)

  ## Initialize empty data frames for each assay (to be filled in the loop)
  cds_ipa_usage <- data.frame()  # IPA usage values
  ipa_reads <- data.frame()      # Raw IPA region read counts
  cds_reads <- data.frame()      # Raw CDS region read counts
  ipa_rpkm <- data.frame()        # IPA region RPKM
  cds_rpkm <- data.frame()        # CDS region RPKM

  ## Loop through each sample and merge their values into the assay matrices
  for (sample in sampleNames) {

    # Read per-sample IPA usage results (CSV)
    x <- read.csv(paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_", sample, "_ipa_usage_atlas.csv"))
    x[is.na(x)] = 0  # Replace NA values with 0 for downstream merging
    print(head(x))   # Print first few rows for progress/debugging

    # Extract and rename columns for each assay, using sample name as column name
    x_cds_ipa_usage <- x[c("X", "cds_ipa_usage")]
    setnames(x_cds_ipa_usage, old = c("cds_ipa_usage"), new = sample)

    x_ipa_reads <- x[c("X", "ipa_reads")]
    setnames(x_ipa_reads, old = c("ipa_reads"), new = sample)

    x_cds_reads <- x[c("X", "cds_reads")]
    setnames(x_cds_reads, old = c("cds_reads"), new = sample)

    x_ipa_rpkm <- x[c("X", "ipa_rpkm")]
    setnames(x_ipa_rpkm, old = c("ipa_rpkm"), new = sample)

    x_cds_rpkm <- x[c("X", "cds_rpkm")]
    setnames(x_cds_rpkm, old = c("cds_rpkm"), new = sample)

    # For the first sample, initialize the assay data frames
    if (nrow(cds_ipa_usage) == 0) {
      cds_ipa_usage <- x_cds_ipa_usage
      ipa_reads     <- x_ipa_reads
      cds_reads     <- x_cds_reads
      ipa_rpkm       <- x_ipa_rpkm
      cds_rpkm       <- x_cds_rpkm
    } else {
      # For subsequent samples, merge by row (IPA event)
      cds_ipa_usage <- merge(cds_ipa_usage, x_cds_ipa_usage)
      ipa_reads     <- merge(ipa_reads, x_ipa_reads)
      cds_reads     <- merge(cds_reads, x_cds_reads)
      ipa_rpkm       <- merge(ipa_rpkm, x_ipa_rpkm)
      cds_rpkm       <- merge(cds_rpkm, x_cds_rpkm)
    }
  }

  # Set row names to IPA event IDs, sort, and remove redundant columns
  rownames(cds_ipa_usage) <- cds_ipa_usage$X
  cds_ipa_usage <- cds_ipa_usage[order(row.names(cds_ipa_usage)), ]
  cds_ipa_usage <- cds_ipa_usage[-1]

  rownames(ipa_reads) <- ipa_reads$X
  ipa_reads <- ipa_reads[order(row.names(ipa_reads)), ]
  ipa_reads <- ipa_reads[-1]

  rownames(cds_reads) <- cds_reads$X
  cds_reads <- cds_reads[order(row.names(cds_reads)), ]
  cds_reads <- cds_reads[-1]

  rownames(ipa_rpkm) <- ipa_rpkm$X
  ipa_rpkm <- ipa_rpkm[order(row.names(ipa_rpkm)), ]
  ipa_rpkm <- ipa_rpkm[-1]

  rownames(cds_rpkm) <- cds_rpkm$X
  cds_rpkm <- cds_rpkm[order(row.names(cds_rpkm)), ]
  cds_rpkm <- cds_rpkm[-1]

  ## Create the rowRanges object for the SummarizedExperiment
  ipa_atlas <- ipa_atlas[order(names(ipa_atlas)), ]
  rowRanges_ipa <- ipa_atlas[names(ipa_atlas) %in% rownames(cds_ipa_usage), ]

  ## Create the colData object (sample metadata)
  colData_ipa <- data.input
  rownames(colData_ipa) <- data.input$NAME

  ## Create the SummarizedExperiment object with all assays and metadata
  ipa.se <- SummarizedExperiment(
    assays = list(
      cds_ipa_usage = cds_ipa_usage,
      ipa_reads = ipa_reads,
      cds_reads = cds_reads,
      ipa_rpkm = ipa_rpkm,
      cds_rpkm = cds_rpkm
    ),
    colData = colData_ipa,
    rowRanges = rowRanges_ipa
  )

  ## Save the SummarizedExperiment object to disk for downstream analysis
  print(paste("Saving SE for : ", sample))
  saveRDS(ipa.se, paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_ipa_usage_se.rds"))

  print("finished")
}
