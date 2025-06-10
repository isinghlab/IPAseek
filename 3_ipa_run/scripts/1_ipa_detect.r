library(GenomicAlignments)
library(GenomicFeatures)
library(data.table)
library(gtools)
library(dplyr)
library(tidyr)

ipa_detect <- function(input.data.path, wd, atlas_name) {


  ## Main function for IPA detection pipeline ##

  # Setup & Initialization

  tryCatch({
    # Create output directory structure
    results.files.dir <- file.path(wd, "pelt", "results", atlas_name)
    slurm.files.dir <- file.path(wd, "pelt", "slurm_submission", atlas_name)
    logs.files.dir <- file.path(wd, "pelt", "logs", atlas_name)
    
    # Create directories with error handling
    safeDirCreate <- function(dir_path) {
      tryCatch({
        if (!dir.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE, showWarnings = TRUE)
        }
      }, error = function(e) {
        stop(paste("Directory creation failed:", dir_path, "\nError:", e$message))
      })
    }
    
    safeDirCreate(results.files.dir)
    safeDirCreate(slurm.files.dir)
    safeDirCreate(logs.files.dir)
    

    # Data Input Handling

    # Read input data with error handling
    data.input <- tryCatch({
      read.delim(input.data.path, sep = "\t", header = TRUE)
    }, error = function(e) {
      stop(paste("Failed to read input data:", input.data.path, "\nError:", e$message))
    })
    

    sampleNames <- as.character(data.input$NAME)
    

    # Sample Processing

    processSample <- function(sample) {

      tryCatch({
        # 3.1 BAM File Handling
        uniq.bam.dir <- as.character(subset(data.input, NAME == sample)$FILE_PATH)
        # safeDirCreate(uniq.bam.dir)
        
        uniq.bam.loc <- file.path(uniq.bam.dir, paste0(sample, "_uniq.bam"))
        reads.dir <- file.path(uniq.bam.dir, paste0(sample, "_reads.rds"))
        
        # 3.2 Filter BAM Files
        if (!file.exists(uniq.bam.loc)) {
          print("No BAM file available")
          print("Use rds object containing Reads")

          if (!file.exists(reads.dir)) { 
            stop("Cannot proceed. No bam or read file available")
            }else{
              # Extract the LIB_SIZE value for the sample
              lib_size_val <- subset(data.input, NAME == sample)$LIB_SIZE

              # Check if LIB_SIZE exists and is not NA or empty
              if (length(lib_size_val) > 0 && !is.na(lib_size_val) && lib_size_val != "") {
                write(lib_size_val, paste0(results.files.dir, "/", sample, "_library_size"))
              }else{
                stop("Cannot proceed. Library Size not provided")
              }
            }

        }else{

          if (!file.exists(reads.dir)){
                    rng_chr<- as(seqinfo(BamFile(uniq.bam.loc)), "GRanges")
                    rng_chr$param<-paste0(seqnames(rng_chr),":1-",end(rng_chr),":",strand(rng_chr))
          
                    bf<-BamFile(uniq.bam.loc)
          
                    #read bam file
                    read_bam<- function(bf) {
                      param <- ScanBamParam(which = rng_chr)
                      reads <- readGAlignments(bf,use.names=TRUE, param = param)
                      }
          
                    print("Read bam")
          
                    reads<-read_bam(bf)
                    saveRDS(reads, reads.dir)}else{
                      reads <- readRDS(reads.dir)
                    }

          # Check if library size has already been computed. If not compute it
          if (!file.exists(paste0(results.files.dir, "/", sample, "_library_size"))) {
          system(paste0("samtools view -c -F 4 ", uniq.bam.loc, ">", results.files.dir, "/", sample, "_library_size"))
          }

        }

        # Expression Data Handling

        exprs.loc <- file.path(
          wd, "2_gene_preprocessing", "3_gene_expression", "results",
          atlas_name, "gene_expression_results", paste0(sample, "_rpkm_df.rds")
        )
        
        exprGenes <- tryCatch({
          readRDS(exprs.loc) %>% 
            as.data.frame() %>% 
            na.omit()
        }, error = function(e) {
          stop(paste("Failed to load expression data for:", sample, "\nError:", e$message))
        })
        

        #Intron Data Processing

        introns.file.dir <- file.path(
          wd, "1_intron_preprocessing", "3_filtering_gobj", 
          "rnhg38_filtered_introns_cds.rds"
        )
        
        rn_introns <- tryCatch({
          readRDS(introns.file.dir) %>% 
            .[names(.) %in% rownames(exprGenes)]
        }, error = function(e) {
          stop(paste("Intron data processing failed for:", sample, "\nError:", e$message))
        })
        
        # Save expressed introns
        exp.files.dir <- file.path(results.files.dir, "exp_introns_results")
        safeDirCreate(exp.files.dir)
        
        tryCatch({
          saveRDS(rn_introns, file.path(exp.files.dir, paste0(sample, ".RDS")))
        }, error = function(e) {
          stop(paste("Failed to save intron data for:", sample, "\nError:", e$message))
        })
        
      

        # Chromosome Processing

        chrs <- unique(seqnames(rn_introns))
        # chrs <- head(chrs)
        
        # Create SLURM submission files
        ipaDetect.loc <- file.path(slurm.files.dir, "ipaDetect.Rdata")
        save(ipaDetect, file = ipaDetect.loc)
        
        # Create sample-specific directories
        sample.dirs <- c(
          file.path(results.files.dir, sample),
          file.path(slurm.files.dir, sample)
        )
        sapply(sample.dirs, safeDirCreate)
        
        # Process each chromosome
        for (chr in chrs) {
          tryCatch({
            # 6.1 Create Genomic Ranges
            GG.locs <- rn_introns[seqnames(rn_introns) == chr]
            GG.locs$output <- file.path(results.files.dir, sample)
            GG.locs$readspath <- reads.dir
            GG.locs$sample <- sample
            GG.locs$wd <- wd
            GG.locs$atlas_name <- atlas_name
            
            # 6.2 Save chromosome data
            chr.data.file <- file.path(
              slurm.files.dir, sample,
              paste0(sample, "_", chr, "GGlocs.Rdata")
            )
            save(GG.locs, file = chr.data.file)
            
            # 6.3 Generate SLURM script
            script.name <- file.path(
              slurm.files.dir, sample,
              paste0(chr, "_ipaDetectRun.R")
            )

            # cat(paste0("\nload (\'",ipaDetect.loc,"')"))

            generateScript <- function() {
              sink(script.name)
              cat(paste0("
                  \nlibrary(GenomicAlignments)
                  \nlibrary(GenomicFeatures)
                  \nlibrary(tidyverse)
                  \nlibrary(dplyr)
                  \nlibrary(data.table)
                  \nload(\'",ipaDetect.loc,"\')
                  \nload(\'",chr.data.file,"\')
                  \nipaDetect(GG.locs)
                  "), fill = TRUE)
              sink()
            }
            
            tryCatch({generateScript()}, error = function(e) {
              stop(paste("Script generation failed for:", chr, "\nError:", e$message))
            })
            
            # 6.4 Submit SLURM job
            bash.file <- file.path(
              slurm.files.dir, sample,
              paste0(chr, 'submit.sh')
            )
            
            slurm.params <- c(
              "#!/bin/bash",
              paste("#SBATCH --job-name=", sample, "_IPA_Detect", sep = ""),
              "#SBATCH --time=01:00:00",
              "#SBATCH --mem=64G",
              "#SBATCH --ntasks=128",
              # "#SBATCH --gres=gpu:a100:1",
              paste("#SBATCH --output=", logs.files.dir, "/", sample, "_IPA_Detect_", chr, sep = ""),
              paste("Rscript --vanilla", script.name)
            )
            
            tryCatch({
              writeLines(slurm.params, bash.file)
              job_output <- system(paste("sbatch", bash.file),  intern = TRUE)

              # Extract only the job ID(s)
              job_id <- regmatches(job_output, regexpr("\\d+", job_output))[[1]]
              write(job_id, "current_jobs.log", append = TRUE)

            }, error = function(e) {
              stop(paste("SLURM submission failed for:", chr, "\nError:", e$message))
            })
            
          }, error = function(e) {
            message(paste("Chromosome processing failed for:", chr, "\nError:", e$message))
          })
        }
        
        message(paste0("Successfully submitted jobs for sample: ", sample))
        
      }, error = function(e) {
        message(paste("Sample processing failed for:", sample, "\nError:", e$message))
        return(NULL)
      })
    }
    
    # Process all samples
    sapply(sampleNames, processSample)
    
  }, error = function(e) {
    message(paste("Critical pipeline error:", e$message))
    return(NULL)
  })
  
  message("Pipeline execution completed")
  invisible(TRUE)
}


################################################################################################################
##########################               Main Function IPA_FINDER                       ########################  
################################################################################################################        

ipaDetect <- function(GG.locs) {
  # Record the start time for performance measurement
  start_time <- Sys.time()
  
  # Helper function for safe file reading
  safeReadRDS <- function(file_path, label = "") {
    tryCatch(
      readRDS(file_path),
      error = function(e) {
        stop(paste("Error reading", label, "file:", file_path, "\n", e$message))
      }
    )
  }
  
  # Helper function for safe directory creation
  safeDirCreate <- function(dir_path) {
    if (!file.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
  
  # Extract parameters
  wd <- GG.locs$wd[1]
  sample_nam <- GG.locs$sample[1]
  result.dir.path <- GG.locs$output[1]
  chr <- seqnames(GG.locs[1])
  
  #  Source external functions with error handling 
  tryCatch({
    source(file.path(wd, "3_ipa_run/scripts/1_ipa_functions.r"))
  }, error = function(e) {
    stop(paste("Error sourcing functions script:", e$message))
  })
  
  #  Coverage calculation 
  coverage.files.dir <- file.path(result.dir.path, "coverage_results")
  coverage.file <- file.path(coverage.files.dir, paste0("cov_", sample_nam, "_", chr, ".RDS"))
  
  if (!file.exists(coverage.file)) {
    message("Coverage file missing. Creating new coverage file...")
    safeDirCreate(coverage.files.dir)
    calCoverage <- tryCatch(
      findCoverage(GG.locs, result.dir.path),
      error = function(e) {
        stop(paste("Error calculating coverage:", e$message))
      }
    )
  } else {
    calCoverage <- safeReadRDS(coverage.file, "coverage")
  }
  
  #  Retention calculation 
  retention.files.dir <- file.path(result.dir.path, "retention_results")
  retention.file <- file.path(retention.files.dir, paste0("retain_introns_", sample_nam, "_", chr, ".RDS"))
  
  if (!file.exists(retention.file)) {
    message("Retention file missing. Creating new retention file...")
    safeDirCreate(retention.files.dir)
    intronRetention <- tryCatch(
      findIntronRetention(calCoverage, sample_nam, retention.files.dir, chr, GG.locs),
      error = function(e) {
        stop(paste("Error calculating intron retention:", e$message))
      }
    )
  } else {
    intronRetention <- safeReadRDS(retention.file, "retention")
  }
  
  #  Execution time 
  exec.time <- Sys.time() - start_time
  message(sprintf("Execution time of coverage and retention calculations: %.2f seconds", as.numeric(exec.time, units="secs")))
  
  #  Return results for further use if needed 
  invisible(list(
    coverage = calCoverage,
    retention = intronRetention,
    execution_time = exec.time
  ))
}




# #################################################################################################################################  
# ###############################  Retreivng intronretention data, Creating sumarized experiment       ############################
# #################################################################################################################################  

# Function to create a table summarizing intron retention information for each sample,
# with robust error handling.
retrieve_intronreten_data <- function(input.data.path, wd, atlas_name) {

  # Define the output directory for intron retention summary results
  intron_retention_data.dir <- file.path(wd, "pelt", "results", atlas_name, "intron_retention_se_results")
  
  # Create the output directory if it does not exist
  if (!dir.exists(intron_retention_data.dir)) {
    dir.create(intron_retention_data.dir, recursive = TRUE)
  }
  
  # Read the input data table, which must contain a column "NAME" for sample names
  data.input <- tryCatch({
    read.delim(input.data.path, sep = "\t", header = TRUE)
  }, error = function(e) {
    stop(paste("Failed to read input data table:", e$message))
  })

  sampleNames <- as.character(data.input$NAME)

  rn_introns <- readRDS(file.path( wd, "1_intron_preprocessing", "3_filtering_gobj", "rnhg38_filtered_introns_cds.rds"))
  expected_chromosomes <- unique(seqnames(rn_introns))
  
  # Process each sample individually with error handling
  sapply(sampleNames, function(sample) {
    tryCatch({
      # Print the current sample name for tracking progress
      print(sample)
      
      # Define the directory containing per-chromosome intron retention files for this sample
      intronret.loc <- file.path(wd, "pelt", "results", atlas_name, sample, "retention_results")
      
      # List all RDS files matching "retain_introns*" (per-chromosome retention results)
      file_list <- list.files(path = intronret.loc, pattern = "retain_introns*", full.names = TRUE)

      present_chromosomes <- gsub(".*retain_introns_.*_(chr[0-9XYM]+)\\.RDS$", "\\1", basename(file_list))

      # Identify missing chromosomes
      missing_chromosomes <- setdiff(expected_chromosomes, present_chromosomes)


      if (length(file_list) == 0) {
        warning(paste("No retention RDS files found for sample:", sample))
        return(NULL)
      } else if (length(missing_chromosomes) > 0) {
        message(paste(
          "Retention files missing for chromosomes:",
          paste(missing_chromosomes, collapse = ", "),
          "in sample:", sample
        ))
      }
      
      # Read all RDS files into a list of data.frames, handling read errors
      data_list <- lapply(file_list, function(f) {
        tryCatch({
          readRDS(f)
        }, error = function(e) {
          warning(paste("Failed to read file:", f, "for sample:", sample, ":", e$message))
          NULL
        })
      })
      
      # Filter out any empty or failed data.frames
      non_empty_data_list <- Filter(function(df) !is.null(df) && nrow(df) > 0, data_list)
      if (length(non_empty_data_list) == 0) {
        warning(paste("All retention files empty or failed for sample:", sample))
        return(NULL)
      }
      
      # Combine all non-empty data.frames into a single data.frame
      mergedat <- tryCatch({
        do.call(rbind, non_empty_data_list)
      }, error = function(e) {
        warning(paste("Failed to merge data for sample:", sample, ":", e$message))
        return(NULL)
      })
      if (is.null(mergedat)) return(NULL)
      
      # Save the merged retention data for this sample
      tryCatch({
        saveRDS(mergedat, file.path(intronret.loc, paste0("retain_all_", sample, ".RDS")))
      }, error = function(e) {
        warning(paste("Failed to save merged retention data for sample:", sample, ":", e$message))
      })
      
      # Create unique intron IDs by combining entrez.id and intron columns
      retained_ids <- tryCatch({
        do.call(paste, c(mergedat[c("entrez.id", "intron")], sep = " "))
      }, error = function(e) {
        warning(paste("Failed to create retained_ids for sample:", sample, ":", e$message))
        return(NULL)
      })
      if (is.null(retained_ids)) return(NULL)
      
      # Load the reference list of all introns from the annotation file
      rn_introns <- tryCatch({
        readRDS(file.path(wd, "1_intron_preprocessing", "3_filtering_gobj", "rnhg38_filtered_introns_cds.rds"))
      }, error = function(e) {
        warning(paste("Failed to load intron annotation for sample:", sample, ":", e$message))
        return(NULL)
      })
      if (is.null(rn_introns)) return(NULL)
      
      # Prepare a data.frame of retained intron IDs and label them as "intron retention"
      retained_ids <- as.data.frame(retained_ids)
      retained_ids$retention <- "intron retention"
      names(retained_ids) <- c("intron_ids", "retention")
      
      # Prepare a data.frame of all possible intron IDs from the annotation
      intron_ids <- as.data.frame(names(rn_introns))
      names(intron_ids) <- "intron_ids"
      
      # Merge the retention status with the full list of introns (left join)
      retention_data <- merge(intron_ids, retained_ids, by = "intron_ids", all.x = TRUE)
      
      # Replace missing values (introns not retained) with "no retention"
      retention_data <- retention_data %>% mutate(across(2, ~ replace_na(.x, "no retention")))
      
      # Rename the retention column to include the sample name for clarity
      colnam <- paste0(sample, "_retention")
      names(retention_data) <- c("intron_ids", colnam)
      
      # Remove duplicate rows, if any
      retention_data <- unique(retention_data)
      
      # Set rownames to intron IDs for easier downstream access
      rownames(retention_data) <- retention_data$intron_ids
      
      # Keep only the retention status column for output
      retention_data <- retention_data[colnam]
      
      # Save the retention data as RDS and CSV for this sample
      tryCatch({
        saveRDS(retention_data, file.path(intron_retention_data.dir, paste0("retention_data_", sample, ".rds")))
        write.csv(retention_data, file.path(intron_retention_data.dir, paste0("retention_data_", sample, ".csv")))
      }, error = function(e) {
        warning(paste("Failed to save retention summary for sample:", sample, ":", e$message))
      })
    }, error = function(e) {
      message(paste("Error processing sample", sample, ":", e$message))
      return(NULL)
    })
  })
}

#################################################################################################################################  
###############################  Function to create table for retention information for all samples  ############################
#################################################################################################################################  

intronret_se <- function(wd, atlas_name) {

  # Define output directory for the combined SummarizedExperiment
  intron_retention.dir <- file.path(wd, "pelt", "results", atlas_name, "results_se")

  # Create the output directory if it does not exist
  if (!dir.exists(intron_retention.dir)) {
    dir.create(intron_retention.dir, recursive = TRUE)
  }

  # Define the path containing per-sample retention data RDS files
  file_path <- file.path(wd, "pelt", "results", atlas_name, "intron_retention_se_results")

  # Try-catch block for the main data loading and processing
  tryCatch({
    # List all RDS files in the retention results directory
    rds_files <- list.files(path = file_path, pattern = "\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) {
      stop(paste("No .rds files found in", file_path))
    }

    # Read all RDS files and combine them column-wise
    retention_list <- lapply(rds_files, function(f) {
      tryCatch({
        readRDS(f)
      }, error = function(e) {
        warning(paste("Failed to read file:", f, ":", e$message))
        NULL
      })
    })
    # Remove any NULL elements from failed reads
    retention_list <- Filter(Negate(is.null), retention_list)
    if (length(retention_list) == 0) {
      stop("No valid retention data could be loaded.")
    }

    # Combine all retention data into a single matrix
    retention <- do.call(cbind, retention_list)

    # Order the rows for consistency using mixedorder (from gtools)
    retention <- retention[mixedorder(rownames(retention)), , drop = FALSE]

    # Create a SummarizedExperiment object from the combined matrix
    se_retention <- SummarizedExperiment(assays = as.matrix(retention))

    # Save the SummarizedExperiment object as an RDS file
    saveRDS(se_retention, file.path(intron_retention.dir, "se_intron_retention.RDS"))

    # Also save the combined matrix as a CSV file
    write.csv(retention, file.path(intron_retention.dir, "se_intron_retention.csv"))

    print("Successfully created and saved combined intron retention SummarizedExperiment.")

  }, error = function(e) {
    message("Error in intronret_se: ", e$message)
  })
}




 
 
 
