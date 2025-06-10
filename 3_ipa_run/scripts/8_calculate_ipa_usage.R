
library(GenomicAlignments)  
library(GenomicFeatures)    
library(dplyr)             
library(tidyr)              
library(data.table)         
library(gtools)             
library(ggplot2)           
library(scales)             

# Function to calculate IPA usage metrics and submit SLURM jobs for parallel processing
calc_ipa_usage <- function(input.data.path, wd, atlas_name) {
  
  # Set up SLURM job submission and log directories
  slurm.files.dir <- file.path(wd, "pelt", "slurm_submission_ipa_usage")
  logs.files.dir <- file.path(wd, "pelt", "logs_ipa_usage")
  
  # Create directories if they don't exist
  dir.create(slurm.files.dir, showWarnings = FALSE)
  dir.create(logs.files.dir, showWarnings = FALSE)

  # Read sample metadata file
  data.input <- read.delim(input.data.path, sep = "\t", header = TRUE)
  sampleNames <- as.character(data.input$NAME)

  # Process each sample in parallel via SLURM
  sapply(sampleNames, function(sample) {
    
    #  Path Configuration 
    filePath <- data.input[data.input$NAME == sample, ]$FILE_PATH
    reads.dir <- file.path(filePath, paste0(sample, "_reads.rds"))
    
    # Validate reads file existence
    if (!dir.exists(filePath)) {
      warning(paste("Missing reads directory for sample:", sample))
    }
    
    # Define IPA atlas and library size paths
    ipa_atlas_path <- paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_full_ipa_atlas_conf.RDS")
    library_size <- read.table(paste0(wd, "/pelt/results/", atlas_name, "/", sample, "_library_size"))$V1

    #  Job Configuration Object 
    paths <- data.frame(
      reads.dir = reads.dir,
      sample = sample,
      condition = data.input[data.input$NAME == sample, "CONDITION"],
      ipa_atlas = ipa_atlas_path,
      wd = wd,
      ipa_atlas_name = atlas_name,
      library_size = library_size
    )

    # SLURM Script Generation 
    # Save function environment
    ipa_usage.loc <- paste0(slurm.files.dir, "/ipa_usage.Rdata")
    save("ipa_usage", file = ipa_usage.loc)
    
    # Create sample-specific directories
    sample.dir <- file.path(slurm.files.dir, sample)
    dir.create(sample.dir, showWarnings = FALSE)
    
    # Generate R script
    script.name <- file.path(sample.dir, "_ipa_usage.R")
    save(paths, file = paste0(sample.dir, "/", sample, "_paths.Rdata"))
    
    sink(file = script.name)
    cat(
      "\nload ('", ipa_usage.loc, "')",
      "\nload ('", sample.dir, "/", sample, "_paths.Rdata')",
      "\nipa_usage(paths)",
      sep = ""
    )
    sink()
    
    #  SLURM Batch Script 
    bash.file.location <- file.path(sample.dir, "submit.sh")
    writeLines(
      c(
        "#!/bin/bash",
        paste("#SBATCH --job-name=", sample, "_ipa_usage", sep = ""),
        "#SBATCH --time=05:00:00",
        "#SBATCH --mem=56G",
        "#SBATCH --ntasks=16",
        paste("#SBATCH --output=", logs.files.dir, "/", sample, "_ipa_usage", sep = ""),
        paste("Rscript --vanilla", script.name)
      ),
      con = bash.file.location
    )
    
    # Submit job
    system(paste("sbatch", bash.file.location))
    print(paste0("Job submitted for sample: ", sample))
  })
}

#################################################################################################################################  
#######################################   Function to calculate IPA usage   #####################################################
#################################################################################################################################  
 
ipa_usage <- function(paths){

    # Extract relevant parameters and file paths from the input 'paths' object
    reads_path      <- paths$reads.dir          # Path to sample's aligned reads (RDS)
    sample          <- paths$sample             # Sample name
    ipa_atlas_path  <- paths$ipa_atlas          # Path to confident IPA atlas (RDS)
    atlas_name      <- paths$ipa_atlas_name     # Atlas name identifier
    library_size    <- paths$library_size       # Library size (for TPM normalization)
    wd              <- paths$wd                 # Working directory

    # Source external script with required functions for UTR3 and CDS calculations
    source(paste0(wd,"/3_ipa_run/scripts/8_calculate_utr3.R"))

    print(paste("calculating IPA usage : ", sample)) # Progress message

    # Load the confident IPA atlas for this sample
    print("Getting IPA_atlas")
    ipa_atlas <- readRDS(ipa_atlas_path)

    # Split the atlas into a list of GRanges objects (one per event)
    ipa_atlas <- split(ipa_atlas, as.factor(ipa_atlas))

    # For each IPA event, calculate CDS information using the external function
    ipa_atlas <- lapply(ipa_atlas, calc_cds, wd)

    # Flatten the list back to a single GRanges object
    ipa_atlas <- unlist(ipa_atlas)
    ipa_atlas <- GRangesList(ipa_atlas)
    ipa_atlas <- unlist(ipa_atlas)

    # Filter out IPA events where CDS start or end is "none"
    ipa_atlas <- ipa_atlas[ipa_atlas$cds_start != "none" | ipa_atlas$cds_end != "none",]

    # Convert the GRanges object to a data frame for easier manipulation
    ipa_atlas_df <- as.data.frame(ipa_atlas)
    # Rename columns for clarity (terminal exon start/end/width)
    setnames(ipa_atlas_df, old = c("start", "end", "width"), new = c("te_start", "te_end", "te_width"))

    print(paste("Making cds atlas"))
    # Create a new GRanges object for the CDS regions using the annotated start/end
    cds_atlas <- makeGRangesFromDataFrame(
        ipa_atlas_df,
        start.field = "cds_start",
        end.field = "cds_end",
        keep.extra.columns = TRUE
    )

    # Load the aligned reads for this sample
    reads <- readRDS(reads_path)

    print(paste("Getting coverage"))

    # Count overlaps (read coverage) for IPA and CDS regions
    ipa_atlas$ipa_reads <- countOverlaps(ipa_atlas, reads, ignore.strand = TRUE)
    ipa_atlas$cds_reads <- countOverlaps(cds_atlas, reads, ignore.strand = TRUE)

    # Calculate widths for IPA and CDS regions
    ipa_atlas$ipa_width <- width(ipa_atlas)
    ipa_atlas$cds_width <- width(cds_atlas)

    print(paste("TPM calculation"))

    # Calculate TPM (Transcripts Per Million) for IPA and CDS regions
    ipa_atlas_df <- as.data.frame(ipa_atlas) %>%
        mutate(
            ipa_tpm = ((ipa_reads / ipa_width) * library_size) / 1000000,
            cds_tpm = ((cds_reads / cds_width) * library_size) / 1000000
        )

    print(paste("IPA usage"))
    # Calculate IPA usage as the fraction of IPA TPM over total (IPA + CDS) TPM
    ipa_atlas_df$cds_ipa_usage <- (ipa_atlas_df$ipa_tpm / (ipa_atlas_df$ipa_tpm + ipa_atlas_df$cds_tpm))

    print(paste("Saving for : ", sample))
    # Save the final annotated data frame as a CSV file for downstream analysis
    write.csv(
        ipa_atlas_df,
        paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_", sample, "_ipa_usage_atlas.csv")
    )

    print("finished") # Indicate completion for this sample

}







