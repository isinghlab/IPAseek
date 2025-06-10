
library(GenomicAlignments)

# Filtering pipeline: prepares and submits per-chromosome filtering jobs for each sample

filtering <- function(input.data.path, wd, atlas_name) {


  # Set up output, slurm, and log directories
  results.files.dir <- file.path(wd, "pelt", "results", atlas_name)
  slurm.files.dir   <- file.path(wd, "pelt", "slurm_submission_filtering", atlas_name)
  logs.files.dir    <- file.path(wd, "pelt", "logs_filtering", atlas_name)

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

  data.input <- data.input[1:100,]

  # Get all sample names
  sampleNames <- as.character(data.input$NAME)

  # Process each sample individually
  sapply(sampleNames, function(sample) {
    tryCatch({
      # Define coverage and retention file locations
      cvg.loc   <- file.path(wd, "pelt", "results", atlas_name)
      reten.loc <- file.path(wd, "pelt", "results", atlas_name, "intron_retention_se_results", paste0("retention_data_", sample, ".rds"))

      # Load retention information
      reten <- tryCatch({
        readRDS(reten.loc)
      }, error = function(e) {
        stop(paste("Failed to read retention file for sample", sample, ":", e$message))
      })
      reten <- as.data.frame(reten)
      reten$id <- rownames(reten)
      nam <- paste0(sample, "_retention")
      reten_id <- reten[which(reten[[nam]] %in% c("no retention")), "id"]

      # Load intron annotation for this sample
      introns.file.dir <- file.path(wd, "pelt", "results", atlas_name, "exp_introns_results", paste0(sample, ".RDS"))
      if (!file.exists(introns.file.dir)) {
        warning(paste("Annotation naming file is missing for sample", sample))
        return(NULL)
      }
      rn_introns <- tryCatch({
        readRDS(introns.file.dir)
      }, error = function(e) {
        stop(paste("Failed to read intron annotation for sample", sample, ":", e$message))
      })
      rn_introns_df <- as.data.frame(rn_introns)
      rn_introns_df$id <- rownames(rn_introns_df)
      rn_introns_df <- rn_introns_df[which(rn_introns_df$id %in% reten_id), ]
      rn_introns <- makeGRangesFromDataFrame(rn_introns_df, keep.extra.columns = TRUE)
      ref <- rn_introns
      chrs <- unique(seqnames(ref))



      # Save the goFilt function for use in the job scripts
      goFilt.loc <- file.path(slurm.files.dir, "goFilt.Rdata")
      tryCatch({
        save("goFilt", file = goFilt.loc)
      }, error = function(e) {
        warning(paste("Failed to save goFilt function for sample", sample, ":", e$message))
      })

      # Ensure sample-specific directories exist
      for (dir in c(file.path(results.files.dir, sample), file.path(slurm.files.dir, sample), logs.files.dir)) {
        if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
      }

      # For each chromosome, prepare and submit a filtering job
      for (chr in chrs) {
        tryCatch({
          # Define coverage file path for this chromosome
          covPath <- file.path(cvg.loc, sample, "coverage_results", paste0("cov_", sample, "_", chr, ".RDS"))
          GG.locs <- ref[which(seqnames(ref) == chr)]
          GG.locs$output <- file.path(results.files.dir, sample)
          GG.locs$sample <- sample
          GG.locs$covPath <- covPath
          GG.locs$wd <- wd
          GG.locs$atlas_name <- atlas_name

          # Save GG.locs object for use in the job script
          results.sample.dir <- file.path(slurm.files.dir, sample)
          GG.locs.rdata.path <- file.path(results.sample.dir, paste0(sample, "_", chr, "GGlocs.Rdata"))
          save(GG.locs, file = GG.locs.rdata.path)

          # Create the R script for this chromosome
          script.name <- file.path(results.sample.dir, paste0(chr, "_goFiltRun.R"))
          sink(file = script.name)
          cat("
                library(GenomicAlignments)
                library(GenomicFeatures)
                library(dplyr)
                library(data.table)
                ")
          cat(paste0("\nload('", goFilt.loc, "')"))
          cat(paste0("\nload('", GG.locs.rdata.path, "')"))
          cat("\ngoFilt(GG.locs)")
          sink()

          # Create the SLURM bash script for this job
          bash.file.location <- file.path(results.sample.dir, paste0(chr, "submit.sh"))
          slurm.jobname <- sprintf("#SBATCH --job-name=%s_IPA_Filtering", sample)
          slurm.time <- "#SBATCH --time=00:30:00"
          slurm.mem <- "#SBATCH --mem=32G"
          slurm.tasks <- "#SBATCH --ntasks=8"
          slurm.output <- paste0("#SBATCH --output=", logs.files.dir, "/", sample, "_IPA_Filtering_", chr)
          sbatch.line <- sprintf("Rscript --vanilla %s", script.name)
          file.conn <- file(bash.file.location)
          writeLines(c("#!/bin/bash", slurm.jobname, slurm.time, slurm.mem, slurm.tasks, slurm.output, sbatch.line), file.conn)
          close(file.conn)

          # Submit the job
          system(paste0("sbatch ", bash.file.location))
        }, error = function(e) {
          warning(paste("Error preparing/submitting job for sample", sample, "chromosome", chr, ":", e$message))
        })
      }

      print(paste0("Running algorithm for sample ", sample, "..."))
      print("Your job is submitted")
    }, error = function(e) {
      warning(paste("Error processing sample", sample, ":", e$message))
    })
  })
}

################################################################################################################
####### Filter introns based on coverage, requiring a continuous stretch of bases above a threshold ############
################################################################################################################        


goFilt <- function(GG.locs){


  # Get the starting time
  start_time <- Sys.time()

  print("Filtering introns")

  # Get the data
  covPath<-GG.locs$covPath[1]
  sample_nam<-GG.locs$sample[1]
  wd<-GG.locs$wd[1]
  atlas_name <- GG.locs$atlas_name[1]

  covbp_dt<-readRDS(paste0(covPath))
  ids<-names(GG.locs)
  covbp_dt<- covbp_dt[ !grepl("cds",covbp_dt$names) , ]
  covbp_dt <- covbp_dt[which(covbp_dt$names %in% ids),]

  rn_introns<-readRDS(file.path(wd,"/pelt/results",atlas_name, "exp_introns_results", paste0(sample_nam,".RDS")))


  #Split the coverage for each intron
  split_list <- split(covbp_dt,covbp_dt$names)
  ls.covg<-lapply(split_list,function(x){
   covls<-x$coverage
  })

  num.reads<- 5
  length.stretch <- 100


  #Filter for the introns having a continuous strecth of 100bp over the given number of reads 

    ls.cov <- lapply(ls.covg, function(x){
  	xx <- ifelse(x >= num.reads, 1, 0)
  	rn <- Rle(xx)
  	df <- data.frame(run = runLength(rn), value = runValue(rn))
  	df.relevant <- subset(df, value == 1)
  	has.coverage <- FALSE
  	if(nrow(df.relevant) > 0){
  		if(any(df.relevant$run >= length.stretch)){
  			has.coverage <- TRUE
  		}
  	}
  	return(has.coverage)
  	})

  tru_introns<-which(unlist(ls.cov))
  tru_introns_lst<-names(tru_introns)
  rn_introns_rel <- rn_introns[which(names(rn_introns) %in% tru_introns_lst)]
  rn_introns_rel <- as.data.frame(rn_introns_rel)


  #Save the results
  path=paste0(GG.locs$output[1],"/filtering_results/")
  dir.create(file.path(path), recursive = TRUE)

  saveRDS(rn_introns_rel,paste0(path, "/filtered_", sample_nam,"_",seqnames(GG.locs[1]), ".RDS"))

  #Calculate the execuation time
  end_time <- Sys.time()
  exec.time<-end_time - start_time
  print(paste0("execution time of filtering:",exec.time))


}