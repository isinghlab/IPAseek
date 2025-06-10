
library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)
library(tidyr)
library(data.table)
library(gtools)
library(changepoint)


# Prepare and submit per-chromosome PELT changepoint jobs for each sample

pelt <- function(input.data.path, wd, atlas_name) {
  # Set up output, slurm, and log directories
  results.files.dir <- file.path(wd, "pelt", "results", atlas_name)
  slurm.files.dir   <- file.path(wd, "pelt", "slurm_submission_pelt", atlas_name)
  logs.files.dir    <- file.path(wd, "pelt", "logs_pelt", atlas_name)

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
  sapply(sampleNames, function(sample) {
    tryCatch({
      # Define coverage, filtering, and retention file locations
      cvg.loc      <- file.path(wd, "pelt", "results", atlas_name)
      filtering.loc<- file.path(wd, "pelt", "results", atlas_name, sample, "filtering_results")
      reten.loc    <- file.path(wd, "pelt", "results", atlas_name, "intron_retention_se_results", paste0("retention_data_", sample, ".rds"))

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


      # Save the runPelt function for use in the job scripts
      runPelt.loc <- file.path(slurm.files.dir, "runPelt.Rdata")
      tryCatch({
        save("runPelt", file = runPelt.loc)
      }, error = function(e) {
        warning(paste("Failed to save runPelt function for sample", sample, ":", e$message))
      })

      # Ensure sample-specific directories exist
      for (dir in c(file.path(results.files.dir, sample), file.path(slurm.files.dir, sample), logs.files.dir)) {
        if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
      }

      # For each chromosome, prepare and submit a PELT job
      for (chr in chrs) {
        tryCatch({
          # Define coverage file path for this chromosome
          covPath <- file.path(cvg.loc, sample, "coverage_results", paste0("cov_", sample, "_", chr, ".RDS"))
          GG.locs <- ref[which(seqnames(ref) == chr)]
          GG.locs$output <- file.path(results.files.dir, sample)
          GG.locs$sample <- sample
          GG.locs$covPath <- covPath
          GG.locs$filteredPath <- filtering.loc
          GG.locs$wd <- wd
          GG.locs$atlas_name <- atlas_name

          # Save GG.locs object for use in the job script
          results.sample.dir <- file.path(slurm.files.dir, sample)
          GG.locs.rdata.path <- file.path(results.sample.dir, paste0(sample, "_", chr, "GGlocs.Rdata"))
          save(GG.locs, file = GG.locs.rdata.path)

          # Create the R script for this chromosome
          script.name <- file.path(results.sample.dir, paste0(chr, "_runPeltRun.R"))
          sink(file = script.name)
          cat("
                library(GenomicAlignments)
                library(GenomicFeatures)
                library(tidyverse)
                library(dplyr)
                library(data.table)
                library(changepoint)
                ")
          cat(paste0("\nload('", runPelt.loc, "')"))
          cat(paste0("\nload('", GG.locs.rdata.path, "')"))
          cat("\nrunPelt(GG.locs)")
          sink()

          # Create the SLURM bash script for this job
          bash.file.location <- file.path(results.sample.dir, paste0(chr, "submit.sh"))
          slurm.jobname <- sprintf("#SBATCH --job-name=%s_runPelt", sample)
          slurm.time <- "#SBATCH --time=9:00:00"
          slurm.mem <- "#SBATCH --mem=8G"
          slurm.tasks <- "#SBATCH --ntasks=8"
          slurm.output <- paste0("#SBATCH --output=", logs.files.dir, "/", sample, "_runPelt_", chr)
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
###################### Detect changepoints in intron coverage using the PELT algorithm #########################
################################################################################################################        


runPelt <- function(GG.locs) {
  # Start timing for performance tracking
  start_time <- Sys.time()
  print("Detecting changepoints")

  # Extract metadata from GG.locs
  covPath    <- GG.locs$covPath[1]
  sample_nam <- GG.locs$sample[1]
  wd         <- GG.locs$wd[1]
  atlas_name <- GG.locs$atlas_name[1]

  # Try to load coverage data
  covbp_dt <- tryCatch({
    readRDS(covPath)
  }, error = function(e) {
    stop(paste("Failed to read coverage RDS:", e$message))
  })
  covbp_dt <- covbp_dt[!grepl("cds", covbp_dt$names), ]

  # Try to load expressed introns for this sample
  rn_introns <- tryCatch({
    readRDS(file.path(wd, "pelt", "results", atlas_name, "exp_introns_results", paste0(sample_nam, ".RDS")))
  }, error = function(e) {
    stop(paste("Failed to read expressed introns RDS:", e$message))
  })

  # Set filtering parameters

    # Set path for storing results and create directory
    path <- paste0(GG.locs$output[1], "/pelt_results/")
    tryCatch({
      dir.create(file.path(path), recursive = TRUE, showWarnings = FALSE)
    }, error = function(e) {
      warning(paste("Failed to create results directory:", e$message))
    })

    # Try to load filtered introns for this chromosome
    filteredPath <- GG.locs$filteredPath[1]
    filteredPath <- paste0(filteredPath,  "/", "filtered_",  sample_nam, "_", seqnames(GG.locs[1]), ".RDS")
    filter_introns <- tryCatch({
      readRDS(filteredPath)
    }, error = function(e) {
      stop(paste("Failed to read filtered introns RDS:", e$message))
    })

    ids <- rownames(filter_introns)
    covbp_dt <- covbp_dt[which(covbp_dt$names %in% ids), ]

    # Split coverage data by intron name
    ls.obj <- split(covbp_dt, covbp_dt$names)

    # Set PELT and filtering parameters
    cost.benefit      <- 100
    intron.length.min <- 1000
    intron.length.max <- 150000
    seg.length        <- 200

    # Filter for introns of appropriate length
    ls.obj.filtered <- ls.obj[sapply(ls.obj, function(x) (nrow(x) >= intron.length.min & nrow(x) <= intron.length.max))]

    # Apply PELT changepoint detection to each intron
    gr <- lapply(ls.obj.filtered, function(dt.gr) {
      tryCatch({
        print(paste("Processing:", dt.gr$names[1], "chromosome:", dt.gr$chr[1]))

        strand.x <- as.character(dt.gr$strand[1])
        coverage <- as.vector(dt.gr$coverage)
        # Remove start/end skews
        coverage <- coverage[4:(length(coverage) - 4)]
        # Reverse direction if antisense strand
        if (strand.x == "-") coverage <- rev(coverage)

        # Run PELT with CROPS penalty range
        cpt.mean.obj <- cpt.mean(
          as.vector(scale(coverage)),
          penalty = "CROPS",
          pen.value = c(100, 10000),
          method = "PELT",
          minseglen = seg.length,
          param.estimates = TRUE
        )

        m.cpts <- cpt.mean.obj@cpts.full
        N <- nrow(m.cpts)
        chk.empty <- which(!is.na(as.vector(m.cpts)))
        vec.penalties <- as.vector(pen.value.full(cpt.mean.obj))

        if (length(vec.penalties) == 2 && N == 1) {
          m.cpts <- rbind(m.cpts, rep(NA, length(m.cpts[1, ])))
        }

        N.new <- nrow(m.cpts)
        Ncpts.real <- 0
        Ncpts.vals <- 0

        if (length(vec.penalties) == N.new && length(chk.empty) > 0) {
          C0 <- vec.penalties[N.new]
          C1 <- vec.penalties[N.new] - vec.penalties[N.new - 1]
          Ncpts.highest.penalty <- length(which(!is.na(m.cpts[N.new, ])))
          diff.pen <- diff(rev(vec.penalties))
          if (Ncpts.highest.penalty > 0) {
            if (C0 > cost.benefit) {
              diff.pen <- c(-C0, diff.pen[1:length(diff.pen)])
              vec.penalties.fraction <- diff.pen[2:length(diff.pen)] / diff.pen[1:(length(diff.pen) - 1)]
              rle.pos <- Rle(ifelse(vec.penalties.fraction >= 0.5, 1, 0))
              rle.l <- runLength(rle.pos)[1]
              rle.v <- runValue(rle.pos)[1]
              if (!is.na(rle.v) && rle.v == 1) {
                Ncpts.vals <- m.cpts[N.new - rle.l, ]
              } else {
                Ncpts.vals <- m.cpts[N.new, ]
              }
              Ncpts.vals <- Ncpts.vals[!is.na(Ncpts.vals)]
              Ncpts.real <- length(Ncpts.vals)
            }
          } else {
            if (C1 > cost.benefit) {
              vec.penalties.fraction <- diff.pen[2:length(diff.pen)] / diff.pen[1:(length(diff.pen) - 1)]
              rle.pos <- Rle(ifelse(vec.penalties.fraction >= 0.5, 1, 0))
              rle.l <- runLength(rle.pos)[1]
              rle.v <- runValue(rle.pos)[1]
              if (!is.na(rle.v) && rle.v == 1) {
                Ncpts.vals <- m.cpts[N.new - rle.l - 1, ]
              } else {
                Ncpts.vals <- m.cpts[N.new - 1, ]
              }
              Ncpts.vals <- Ncpts.vals[!is.na(Ncpts.vals)]
              Ncpts.real <- length(Ncpts.vals)
            }
          }
        }

        # Map change points to genomic coordinates
        if (strand.x == "+") {
          Ncpts.chr.vals <- dt.gr[1, ]$start + Ncpts.vals
        } else {
          Ncpts.chr.vals <- dt.gr[nrow(dt.gr), ]$start - Ncpts.vals
        }
        end <- Ncpts.chr.vals

        gr.chng <- GRanges(
          seqnames = as.character(dt.gr$chr[1L]),
          ranges = IRanges(start = Ncpts.chr.vals, end = end),
          strand = strand.x,
          entrez.id = as.character(dt.gr$gene[1L]),
          id = as.character(dt.gr$names[1L]),
          Nchng.pts = Ncpts.real
        )
        return(gr.chng)
      }, error = function(e) {
        warning(paste("Error in PELT for intron", dt.gr$names[1], ":", e$message))
        return(NULL)
      })
    })

    # Combine results and save
    gr.new <- do.call(c, unlist(gr, use.names = FALSE))
    gr.new2 <- as.data.frame(gr.new)

    tryCatch({
      saveRDS(gr.new, paste0(path, "/ipa_final_", sample_nam, "_filtered_",  seqnames(GG.locs[1]), ".RDS"))
      write.csv(gr.new, paste0(path, "/ipa_final_", sample_nam, "_filtered_",  seqnames(GG.locs[1]), ".csv"))
      saveRDS(gr.new2, paste0(path, "/ipa_final_df", sample_nam, "_filtered_",  seqnames(GG.locs[1]), ".rds"))
      print("Done: change points detected")
    }, error = function(e) {
      warning(paste("Failed to save changepoint results:", e$message))
    })

    # Calculate execution time for this threshold
    end_time <- Sys.time()
    exec.time <- end_time - start_time
    print(paste0("execution time of change point analysis: ", exec.time))

}
