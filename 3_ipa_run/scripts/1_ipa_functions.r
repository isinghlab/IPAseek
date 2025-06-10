################################################################################################################  
########################## findCoverage (Calculate genomic coverage over each base pair)########################
################################################################################################################  

findCoverage <- function(GG.locs, result.dir.path) {

  coverage.files.dir <- file.path(result.dir.path, "coverage_results")

  # Start the system time to calculate the processing time
  start_time <- Sys.time()

  print("calculating genome coverage")

  # Extract working directory from input object
  wd <- GG.locs$wd[1]

  # Load the genome annotation file with error handling
  rdsSource <- tryCatch({
    readRDS(paste0(wd, "/1_intron_preprocessing/1_flatten_genome/hg38_annotated_numbered_cds.rds"))
  }, error = function(e) {
    stop(paste("Failed to read genome annotation:", e$message))
  })

  # Get the path for the BAM file
  reads_path <- GG.locs$readspath[1]

  # Get the sample name
  sample_nam <- GG.locs$sample[1]

  # Get the chromosome number and convert to proper format
  seqnames <- tryCatch({
    as.data.table(seqnames(GG.locs[1]))
  }, error = function(e) {
    stop(paste("Failed to extract chromosome name:", e$message))
  })
  seqnames <- paste("chr", toString(seqnames), sep = "")
  if (seqnames == "chr23") seqnames <- "chrX"
  if (seqnames == "chr24") seqnames <- "chrY"

  # # Get the reads for the chromosome
  reads <- readRDS(reads_path)


  # Get introns in the expressed genes (from GG.locs)
  rn_introns <- GG.locs

  # Find genome annotations overlapping end-to-end with the filtered introns
  filt_reg <- tryCatch({
    findOverlaps(GG.locs, rdsSource, type = "end")
  }, error = function(e) {
    stop(paste("Failed to find overlaps:", e$message))
  })
  index <- subjectHits(filt_reg)
  strand <- as.character(strand(rdsSource[index]))

  # Get the index of genome annotation before and after the intron
  index_up <- ifelse(strand == "+", index + 1, index - 1)
  index_dwn <- ifelse(strand == "+", index - 1, index + 1)

  # Combine the indexes and get the elements from the genome annotations
  sel_index <- c(rbind(index_dwn, index, index_up))
  selc <- rdsSource[sel_index]

  # Assign the names for up/downstream and current regions
  nam <- names(rn_introns)
  nam_cds_up <- paste(nam, "_cds_up", sep = '')
  nam_cds_dwn <- paste(nam, "_cds_down", sep = '')
  nams <- c(rbind(nam_cds_up, nam, nam_cds_dwn))
  names(selc) <- nams

  # Group the introns by entrez id and intron number
  grp <- rep(nam, each = 3)
  selc$group <- grp

  # Divide the region into single base pairs (tiles)
  tiles <- tile(selc, width = 1)
  tiles2 <- unlist(unname(tiles))

  # Get coverage over each base
  hits.df <- as.data.frame(countOverlaps(tiles2, reads, ignore.strand = TRUE))
  covbp_dt <- as.data.table(hits.df)
  setnames(covbp_dt, "countOverlaps(tiles2, reads, ignore.strand = TRUE)", "coverage")

  # Create the coverage data table with annotation columns
  rng <- ranges(tiles2)
  covbp_dt$names <- names(rng)
  chr <- as.data.table(seqnames(tiles2[1]))
  st_pos <- start(rng)

  sp_nam <- strsplit(names(rng), ' ')
  gene <- sapply(sp_nam, function(x) x[1])
  intron <- sapply(sp_nam, function(x) x[2])
  covbp_dt$start <- st_pos
  covbp_dt$names <- names(rng)
  covbp_dt$gene <- gene
  covbp_dt$intron <- intron
  strand_dt <- as.data.table(strand(tiles2))
  setnames(strand_dt, "strand")
  chrom <- as.data.table(seqnames(tiles2))
  setnames(chrom, "chr")
  covbp_dt <- cbind(covbp_dt, strand_dt)
  covbp_dt <- cbind(covbp_dt, chrom)

  # Try to save the coverage data table as an RDS file
  tryCatch({
    saveRDS(covbp_dt, paste0(coverage.files.dir, "/cov_", sample_nam, "_", seqnames, ".RDS"))
  }, error = function(e) {
    warning(paste("Failed to save coverage RDS file:", e$message))
  })

  # Print completion message and execution time
  print("Done: genome coverage calculated")
  end_time <- Sys.time()
  exec.time <- end_time - start_time
  print(paste0("execution time of coverage calculation:", exec.time))

  # Return the coverage data table
  return(covbp_dt)
}

 
#################################################################################################################################  
##########################   intronRetention(Identify retained introns and filter them)   #######################################
#################################################################################################################################  


findIntronRetention <- function(covbp_dt, sample_nam, retention.files.dir, chr, GG.locs) {
  # Start the system time to calculate the processing time
  start_time <- Sys.time()

  # Extract working directory from input object
  wd <- GG.locs$wd[1]

  print("Detecting intron retention")
  tryCatch({
    # Extract atlas name
    atlas_name <- GG.locs$atlas_name[1]

    # Fetch the file containing the filtered introns
    filt_genes <- tryCatch({
      readRDS(file.path(wd, "pelt", "results", atlas_name, "exp_introns_results", paste0(sample_nam, ".RDS")))
    }, error = function(e) {
      stop(paste("Failed to read filtered introns:", e$message))
    })

    # Split the names column to create a group column
    sp_nam <- strsplit(covbp_dt$names, '_')
    grp <- sapply(sp_nam, function(x) x[1])
    covbp_dt$group <- grp

    # Extract output directory
    results.files.dir <- GG.locs$output[1]
    # Extract location of uniquely mapped reads bam file
    # uniq.bam.loc <- GG.locs$bampath[1]

  
    # Read library size.
    library_size <- tryCatch({
      read.table(paste0(wd, "/pelt/results/", atlas_name, "/", sample_nam, "_library_size"))$V1
    }, error = function(e) {
      stop(paste("Failed to read library size file:", e$message))
    })

    # Calculate TPM. The tail and head are used to find the start and end pos of the features which are then used to calculate length
    covbp_dt <- covbp_dt %>% group_by(names) %>% mutate(tpm = ((coverage / (tail(start, 1) - head(start, 1))) * library_size) / 1000000)

    ####################### Filter 1: check if 3 reads spanning upstream downstream exon intron junction #############################
    print("filter1")

    # junction_count_up<-covbp_dt %>% group_by(names) %>% summarize(junction_count_up=head(coverage,1))
    # junction_count_down<-covbp_dt %>% group_by(names) %>% summarize(junction_count_down=tail(coverage,1))

    # Get the number of junction reads. The strand argument is used to perform different ops depending on which strand the feature is on
    junction_count_up <- covbp_dt %>% group_by(names) %>% reframe(junction_count_up = ifelse(as.character(strand) == "+", tail(coverage, 1), head(coverage, 1)))
    junction_count_down <- covbp_dt %>% group_by(names) %>% reframe(junction_count_down = ifelse(as.character(strand) == "+", head(coverage, 1), tail(coverage, 1)))

    # Use the not grepl function to remove rows where up is not found in the junction data
    dwn_cds_junction <- unique(junction_count_up[!grepl("down", junction_count_up$names), ])
    up_cds_junction <- unique(junction_count_down[!grepl("up", junction_count_down$names), ])

    # dwn_cds_junction<-junction_count_down[!grepl("up",junction_count_up$names) , ]
    # up_cds_junction<-junction_count_up[!grepl("down",junction_count_down$names) , ]

    # Add a column group which groups the features from up or downstream into groups
    sp_nam <- strsplit(dwn_cds_junction$names, '_')
    grp <- sapply(sp_nam, function(x) x[1])
    dwn_cds_junction$group <- grp

    # Add a column group which groups the features from up or downstream into groups
    sp_nam <- strsplit(up_cds_junction$names, '_')
    grp <- sapply(sp_nam, function(x) x[1])
    up_cds_junction$group <- grp

    # If junction count is less than equal to three then no retention downstream or up, else intron retention.
    dwn_span <- dwn_cds_junction %>% group_by(group) %>% summarise(retention = if (any(junction_count_up <= 3)) {
      "no retention_dwn"
    } else {
      "intron retention_dwn"
    })

    # If junction count is less than equal to three then no retention upstream or up, else intron retention.
    up_span <- up_cds_junction %>% group_by(group) %>% summarise(retention = if (any(junction_count_down <= 3)) {
      "no retention_up"
    } else {
      "intron retention_up"
    })

    # Rbind up and downstream spanning reads and create a logical to find introns
    filter1 <- rbind(dwn_span, up_span)
    filter1_introns <- filter1[grepl("intron", filter1$retention), ]

    # Get the names of introns which are only retained by filter1
    filter1_introns_names <- filter1_introns$group
    passed_filter1_introns <- covbp_dt[covbp_dt$group %in% filter1_introns_names, ]

    # Save the features and coverage of features retained by first filter
    saveRDS(passed_filter1_introns, paste0(retention.files.dir, "/IR_filter1_", sample_nam, "_", chr, ".RDS"))

    ####################### Filter 2: Check at least 50% of the intron length should be covered by more than 1TPM unique reads ################
    print("filter2")

    # Select for the regions in cov bp which are not upstream or downstream. These are the intron regions for which we calculate TPMs
    passed_filter1_introns2 <- passed_filter1_introns[(!grepl("down", passed_filter1_introns$names) & !grepl("up", passed_filter1_introns$names)), ]

    # Calculate the total number of reads, number of reads more than 1, and number of reads more than 0
    unique_reads <- passed_filter1_introns2 %>% group_by(names) %>% summarise(total_count = n(), read_count = sum(tpm > 1), read_count_0 = sum(tpm > 0))
    # Calculate the percentage of reads in each group which are more than zero or more than 1.
    unique_reads <- unique_reads %>% group_by(names) %>% mutate(read_percentage = (read_count / total_count) * 100, read_percentage_0 = (read_count_0 / total_count) * 100)
    # For features which are more than 85 percent we say that the intron is retained else not.
    filter2 <- unique_reads %>% group_by(names) %>% summarise(retention = if (any(read_percentage_0 >= 85)) {
      "intron retention"
    } else {
      "no retention"
    })

    # Subset for the introns whose names match the intron names from filter 2
    filter2_introns <- filter2[grepl("intron", filter2$retention), ]

    # Create group IDs based on filter2 output
    sp_nam <- strsplit(filter2_introns$names, '_')
    grp <- sapply(sp_nam, function(x) x[1])
    filter2_introns$group <- grp

    # From the above extracted groupIDs get the names for group
    filter2_introns_names <- filter2_introns$group
    # Select for the group names. Which are in filter 2 introns.
    passed_filter2_introns <- covbp_dt[covbp_dt$group %in% filter2_introns_names, ]

    # Save the retained intron files.
    saveRDS(passed_filter2_introns, paste0(retention.files.dir, "/IR_filter2_", sample_nam, "_", chr, ".RDS"))

    ################################ Filter 3:Check if median coverage over the flanking exons is more than 1TPM  #######################
    print("filter3")

    # The following two extracts the regions in bp cov which up or downstream are present and puts them into their respective tables.
    flanking_exons_up <- passed_filter2_introns[grepl("up", passed_filter2_introns$names), ]
    flanking_exons_dwn <- passed_filter2_introns[grepl("down", passed_filter2_introns$names), ]

    # Rbind the up and downstream data
    flanking_exons <- rbind(flanking_exons_up, flanking_exons_dwn)

    # Convert the table to data frame summarize it to get the median coverage, TPM, mean coverage and mean TPM for each bp in CDS
    cov_flanking_exons <- as.data.frame(flanking_exons) %>% group_by(names) %>% summarise(medCVG = median(coverage), medTPM = median(tpm), meanCVG = mean(coverage), meanTPM = mean(tpm))

    # Create a grouping column again based on the name
    sp_nam <- strsplit(cov_flanking_exons$names, '_')
    grp <- sapply(sp_nam, function(x) x[1])
    cov_flanking_exons$group <- grp

    # Calculate using any if medTPM is less than zero otherwise the intron is retained
    covbp_flanking_exons <- cov_flanking_exons %>% group_by(group) %>% summarise(retention = if (any(medTPM <= 1)) {
      "no retention"
    } else {
      "intron retention"
    })

    # The the intron names are are retained
    filter3_introns_names <- covbp_flanking_exons$group
    # Select for the group for which intron names are what was extracted in above line.
    passed_filter3_introns <- covbp_dt[covbp_dt$group %in% filter3_introns_names, ]

    # Save all the files which pass this TPM calculation
    saveRDS(passed_filter3_introns, paste0(retention.files.dir, "/IR_filter3_", sample_nam, "_", chr, ".RDS"))

    ####### Filter 4:Check if the ratio of median coverage over the intron to median coverage of the upstream exon was at least 10%  #######
    print("filter4")

    # Calculate the median value for those introns that are passed from above line
    medcvg <- passed_filter3_introns %>% group_by(names) %>% summarise(medCVG = median(coverage))

    # CReate groupID column
    sp_nam <- strsplit(medcvg$names, '_')
    grp <- sapply(sp_nam, function(x) x[1])
    medcvg$group <- grp

    # Get the retention ratios
    medcvg2 <- medcvg %>% group_by(group) %>% mutate(ratio = abs(first(medCVG) / medCVG))
    # This selects the CDS regions only to find median
    medcvg3 <- medcvg2[grepl("cds", medcvg2$names), ]

    # Set infinte values to 0
    medcvg3$ratio[which(!is.finite(medcvg3$ratio))] <- 0

    # Retain the regions for which any is more than equal to 10 or more is retained other wise not retained.
    filter4 <- medcvg3 %>% group_by(group) %>% summarise(retention = if (any(ratio <= 0.10)) {
      "no retention"
    } else {
      "intron retention"
    })

    # Svae data
    saveRDS(filter4, paste0(retention.files.dir, "/IR_filter4_", sample_nam, "_", chr, ".RDS"))

    ###############################  Selecting retained introns and saving the retention ratios in file  #################################

    # Select group of introns whcih are actually retained
    filter4_introns <- filter4[grepl("intron retention", filter4$retention), ]
    ids_retintron <- filter4_introns$group

    # Get the retained introns based on filrered genes and group ids
    retained_introns <- filt_genes[which(names(filt_genes) %in% ids_retintron)]

    # From median table get those whihc were passed or in those ids.
    retained_introns_ratios <- medcvg2[medcvg2$group %in% ids_retintron, ]
    # assign the values in data frame
    retained_intron <- as.data.frame(retained_introns)
    # Slect the median values fo those groupIDs and rename the columns
    medcvg_intron <- retained_introns_ratios[!grepl("cds", retained_introns_ratios$names), ]
    colnames(medcvg_intron) <- paste(colnames(medcvg_intron), "intron", sep = "_")
    medcvg_intron <- as.data.frame(medcvg_intron)
    # Get the downstream values of this as well and relabel it
    medcvg_down <- retained_introns_ratios[grepl("down", retained_introns_ratios$names), ]
    colnames(medcvg_down) <- paste(colnames(medcvg_down), "down", sep = "_")
    medcvg_down <- as.data.frame(medcvg_down)
    # Finally get the upstream values and get it in one place
    medcvg_up <- retained_introns_ratios[grepl("up", retained_introns_ratios$names), ]
    colnames(medcvg_up) <- paste(colnames(medcvg_up), "up", sep = "_")
    medcvg_up <- as.data.frame(medcvg_up)
    # Now cbind to form everything in a single matrix
    retained_introns_data <- cbind(retained_intron, medcvg_intron, medcvg_up, medcvg_down)
    retained_introns_data <- as.data.frame(retained_introns_data)
    retained_introns$medCVG_intron <- retained_introns_data$medCVG_intron
    retained_introns$medCVG_up <- retained_introns_data$medCVG_up
    retained_introns$ratio_up <- retained_introns_data$ratio_up
    retained_introns$medCVG_down <- retained_introns_data$medCVG_down
    retained_introns$ratio_down <- retained_introns_data$ratio_down
    retained_introns <- as.data.frame(retained_introns)
    print("saving")

    # Save the retained reads here
    saveRDS(retained_introns, paste0(retention.files.dir, "/retain_introns_", sample_nam, "_", chr, ".RDS"))

    retained_introns <- data.frame(lapply(retained_introns, as.character), stringsAsFactors = FALSE)
    retained_introns_data <- as.data.frame(retained_introns_data)
    retained_introns_data <- data.frame(lapply(retained_introns_data, as.character), stringsAsFactors = FALSE)
    write.csv(retained_introns, paste0(retention.files.dir, "/retained_introns_sel_", sample_nam, "_", chr, ".csv"))
    write.csv(retained_introns_data, paste0(retention.files.dir, "/retained_introns_data_", sample_nam, "_", chr, ".csv"))
    return(retained_introns)
  }, error = function(e) {
    stop(paste("Error in findIntronRetention function: ", e$message))
  })

  # Compute time
  print("Done: intron retention detected")
  end_time <- Sys.time()
  exec.time <- end_time - start_time
  print(paste0("execution time of intron retention detection:", exec.time))
}




 
 
 
