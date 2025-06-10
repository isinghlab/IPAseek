
library(GenomicAlignments)
library(dplyr)             
library(tidyverse)         
library(data.table)        

#################################################################################################################################  
################### function to merge changepoint (cpt) results across multiple samples   #######################################
#################################################################################################################################  

merge_cpts <- function(input.data.path, wd, atlas_name) {

   # Read the input metadata file (tab-delimited), which contains sample names and BAM file locations
   data.input <- read.delim(input.data.path, sep="\t", header=TRUE)

   # Extract sample names as a character vector
   sampleNames <- as.character(data.input$NAME)

   # Iterate over each sample and perform merging
   sapply(sampleNames, function(sample_nam) {

      # Construct the path to the directory containing changepoint results for this sample
      cpt.loc <- paste0(wd, "/pelt/results/", atlas_name, "/", sample_nam, "/filter_cpts_results/")
      # Define the output file prefix for this sample
      out <- paste0(cpt.loc, "/", sample_nam, "_")

      # Merge all filtered changepoint expression CSVs matching "filter_cpts_tpm_exprs*"
      # These files contain expression-filtered changepoint events (TPM-filtered)
      cpt1 <- do.call(rbind, lapply(
        list.files(path = cpt.loc, pattern = "filter_cpts_tpm_exprs*", full.names = TRUE),
        read.csv
      ))

      # Merge all changepoint CSVs matching "filter_cpts_all_"
      # These files contain all detected changepoint events for the sample
      cpt_all <- do.call(rbind, lapply(
        list.files(path = cpt.loc, pattern = "filter_cpts_all_", full.names = TRUE),
        read.csv
      ))
      # Remove duplicate rows to ensure uniqueness
      cpt_all <- unique(cpt_all)

      # Prepare final expression-filtered changepoint table, removing duplicates
      cpt_exprs <- cpt1
      cpt_exprs <- unique(cpt_exprs)

      # Write the merged, unique changepoint tables to CSV files for downstream analysis
      write.csv(cpt_all, paste0(out, "cpt_all.csv"))
      write.csv(cpt_exprs, paste0(out, "cpt_exprs_all.csv"))

      # Print progress message for this sample
      print("finished merging")
   })
}


