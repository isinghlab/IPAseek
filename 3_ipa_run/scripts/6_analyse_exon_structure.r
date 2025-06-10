

library(GenomicAlignments) 
library(GenomicFeatures)    
library(dplyr)              
library(tidyr)              
library(data.table)         
library(gtools)             
library(ggplot2)            
library(scales)            

# Main processing function 
filter_te_run <- function(input.data.path, wd, atlas_name) {


  # Set up output, slurm, and log directories
  results.files.dir <- file.path(wd, "pelt", "results", atlas_name)
  slurm.files.dir   <- file.path(wd, "pelt", "slurm_submission_filt_te", atlas_name)
  logs.files.dir    <- file.path(wd, "pelt", "logs_filt_te", atlas_name)

  # Create directories if they do not exist
  for (dir in c(results.files.dir, slurm.files.dir, logs.files.dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }

  # Read the sample metadata table with error handling
  data.input <- tryCatch({
    read.delim(input.data.path, sep = "\t", header = TRUE)
  }, error = function(e) {
    stop(paste("Failed to read input data table:", e$message))
  })

  # Get all sample names
  sampleNames <- as.character(data.input$NAME)

  # Process each sample individually
  sapply(sampleNames, function(sample_nam) {

    tryCatch({
        # Get path for unique reads file
        reads_path <- as.character(subset(data.input, NAME == sample_nam)$FILE_PATH)

        reads.dir <- file.path(reads_path, paste0(sample_nam, "_reads.rds"))

        # If uniq reads directory does not exist, create it and print a message
        if (!dir.exists(reads_path)) {
          warning(paste("Reads file not found."))
          }

      # Create a file for library size if it doesn't exist
      sample.results.dir <- file.path(results.files.dir, sample_nam)
      if (!dir.exists(sample.results.dir)) {
        dir.create(sample.results.dir, recursive = TRUE)
      }

      # If library size file doesn't exist, calculate and store it
      libsize.file <- file.path(results.files.dir, paste0(sample_nam, "_library_size"))


      # Read library size with error handling
      library_size <- tryCatch({
        read.table(libsize.file)$V1
      }, error = function(e) {
        stop(paste("Failed to read library size file for sample", sample_nam, ":", e$message))
      })

      # Get change point and coverage locations
      cpts.path <- file.path(wd, "pelt", "results", atlas_name, sample_nam, "filter_cpts_results")
      cvg.loc   <- file.path(wd, "pelt", "results", atlas_name, sample_nam, "coverage_results")

      # Get annotation file location
      anno.loc <- file.path(wd, "1_intron_preprocessing", "3_filtering_gobj", "rnhg38_filtered_introns_cds.rds")
      if (!file.exists(anno.loc)) {
        warning("Annotation file is missing!")
        return(NULL)
      }

      # Read annotation file with error handling
      rdsSource <- tryCatch({
        readRDS(anno.loc)
      }, error = function(e) {
        stop(paste("Failed to read annotation file:", e$message))
      })

      # Retrieve the annotations and intron data
      introns.file.dir <- file.path(wd, "pelt", "results", atlas_name, "exp_introns_results", paste0(sample_nam, ".RDS"))
      if (!file.exists(introns.file.dir)) {
        warning("Annotation naming file is missing!")
        return(NULL)
      }

      rn_introns <- tryCatch({
        readRDS(introns.file.dir)
      }, error = function(e) {
        stop(paste("Failed to read intron annotation for sample", sample_nam, ":", e$message))
      })
      ref <- rn_introns

      # Save the filter_te function for use in the job scripts
      filter_te.loc <- file.path(slurm.files.dir, "filter_te.Rdata")
      tryCatch({
        save("filter_te", file = filter_te.loc)
      }, error = function(e) {
        warning(paste("Failed to save filter_te function for sample", sample_nam, ":", e$message))
      })

      # Ensure sample-specific directories exist
      for (dir in c(file.path(slurm.files.dir, sample_nam), logs.files.dir)) {
        if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
      }

      # Prepare GG.locs object with all relevant info
      GG.locs <- ref
      GG.locs$output <- sample.results.dir
      GG.locs$readspath <- reads.dir
      GG.locs$cptspath <- cpts.path
      GG.locs$sample_nam <- sample_nam
      GG.locs$cvgPath <- cvg.loc
      GG.locs$library_size <- library_size
      GG.locs$wd <- wd

      # Save GG.locs object for use in the job script
      results.sample.dir <- file.path(slurm.files.dir, sample_nam)
      GG.locs.rdata.path <- file.path(results.sample.dir, paste0(sample_nam, "_", sample_nam, "GGlocs.Rdata"))
      save(GG.locs, file = GG.locs.rdata.path)

      # Create the R script for this sample
      script.name <- file.path(results.sample.dir, paste0(sample_nam, "_filter_teRun.R"))
      sink(file = script.name)
      cat("
            library(GenomicAlignments)
            library(GenomicFeatures)
            library(tidyverse)
            library(dplyr)
            library(DESeq2)
            library(data.table)
            ")
      cat(paste0("\nload('", filter_te.loc, "')"))
      cat(paste0("\nload('", GG.locs.rdata.path, "')"))
      cat("\nfilter_te(GG.locs)")
      sink()

      # Create the SLURM bash script for this job
      bash.file.location <- file.path(results.sample.dir, paste0(sample_nam, "submit.sh"))
      slurm.jobname <- sprintf("#SBATCH --job-name=%s_filter_te_", sample_nam)
      slurm.time <- "#SBATCH --time=01:30:00"
      slurm.mem <- "#SBATCH --mem=56G"
      slurm.tasks <- "#SBATCH --ntasks=8"
      slurm.output <- paste0("#SBATCH --output=", logs.files.dir, "/", sample_nam, "_filter_te_", sample_nam)
      sbatch.line <- sprintf("Rscript --vanilla %s", script.name)
      file.conn <- file(bash.file.location)
      writeLines(c("#!/bin/bash", slurm.jobname, slurm.time, slurm.mem, slurm.tasks, slurm.output, sbatch.line), file.conn)
      close(file.conn)

      # Submit the job
      system(paste0("sbatch ", bash.file.location))

      print(paste0("Running algorithm for sample ", sample_nam, "..."))
      print("Your job is submitted")
    }, error = function(e) {
      warning(paste("Error processing sample", sample_nam, ":", e$message))
    })
  })
}

#################################################################################################################################  
#################################   Quantify terminal exon (TE) expression)   ###################################################
#################################################################################################################################  

filter_te <- function(GG.locs) {

  #  Extract sample-specific paths and parameters 
  reads_path    <- GG.locs$readspath[1]         # Path to reads file
  cpts_path   <- GG.locs$cptspath[1]        # Path to changepoint results directory
  library_size <- GG.locs$library_size[1]   # Library size for normalization
  sample      <- GG.locs$sample_nam[1]      # Sample name/ID
  wd          <- GG.locs$wd[1]              # Working directory
  cvg_path    <- GG.locs$cvgPath[1]         # Coverage path (not used here)
  outPath     <- GG.locs$output[1]          # Output directory

  print(sample) # Log current sample being processed

  #  Prepare output directory for results 
  path <- file.path(outPath, "exon_exprs_results")
  ipa.dir <- path
  if (!dir.exists(ipa.dir)) {
    dir.create(ipa.dir, recursive = TRUE)
  }

  #  Load changepoint (IPA candidate) data for this sample 
  cpts.loc <- file.path(cpts_path,  paste0(sample, "_cpt_exprs_all.csv"))
  cpts.file <- read.csv(cpts.loc)

  #  Adjust exon boundaries for composite IPA events 
  # For '+' strand: set exon.start to cds.end if appropriate
  # For '-' strand: set exon.end to cds.start if appropriate
  cpts.file_new <- cpts.file %>%
    mutate(
      exon.start = ifelse(
        cds.end != 0 & strand == "+" & ipa_sel == "composite_ipa" &
          exon.start < cds.start & exon.end > cds.end,
        cds.end, exon.start
      ),
      exon.end = ifelse(
        cds.start != 0 & strand == "-" & ipa_sel == "composite_ipa" &
          exon.start < cds.start & exon.end > cds.end,
        cds.start, exon.end
      )
    )

  #  Convert data frame to GRanges object for overlap analysis 
  exon.obj <- makeGRangesFromDataFrame(
    cpts.file_new,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field = "exon.start",
    end.field = "exon.end"
  )
  exon.obj$names <- paste0(exon.obj$id, "_", exon.obj$ipa_sel) # Unique names

  print("Read bam") # Log BAM reading step

  reads <- readRDS(reads_path) # Read alignments

  #  Quantify coverage (number of reads overlapping each exon region) 
  exon.obj$te_coverage <- countOverlaps(exon.obj, reads)

  #  Convert GRanges to data.table for further manipulation 
  exon.obj.df <- as.data.table(exon.obj)

  #  Calculate TPM (Transcripts Per Million) for each region 
  # TPM = (read count / region width) * (library size / 1,000,000)
  exon.obj.df <- as.data.frame(exon.obj.df) %>%
    mutate(te_tpm = ((te_coverage / abs(end - start)) * library_size) / 1e6)

  # Add sample name as a source column
  exon.obj.df$source <- sample

  #  Filter out low-expression regions (TPM < 1) 
  exon.obj.df <- exon.obj.df[exon.obj.df$te_tpm >= 1, ]

  #  Save results as CSV for downstream analysis 
  write.csv(
    exon.obj.df,
    file = file.path(path, paste0("te_expression", sample, ".csv")),
    row.names = FALSE
  )

  print("finished") # Log completion for this sample
}
