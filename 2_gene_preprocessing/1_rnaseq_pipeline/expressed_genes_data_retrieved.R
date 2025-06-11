library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)
library(tidyr)
library(data.table)
library(gtools)

retrieve_geneexpr_data <- function(input.data.path, wd, sample.name) {
  # Function to retrieve and process gene expression data for a set of samples

  # Define the directory to store gene expression results
  gene_expression_data.dir <- file.path(
    paste0(wd, "/2_gene_preprocessing/3_gene_expression/results/", sample.name, "/gene_expression_results")
  )
  
  # Create the gene expression results directory if it doesn't exist
  if (!dir.exists(gene_expression_data.dir)) {
    dir.create(gene_expression_data.dir)
  }
  
  # Define the directory to store summarized experiment results
  gene_expression_results.dir <- file.path(
    paste0(wd, "/2_gene_preprocessing/3_gene_expression/results/", sample.name, "/gene_expression_se_results")
  )
  
  # Create the summarized experiment results directory if it doesn't exist
  if (!dir.exists(gene_expression_results.dir)) {
    dir.create(gene_expression_results.dir)
  }
  
  # Read the input data table containing sample names and BAM file locations
  data.input <- read.delim(input.data.path, sep = "\t", header = TRUE)
  
  # Extract the list of sample names from the input table
  sampleNames <- as.character(data.input$NAME)
  
  # Process each sample individually
  sapply(sampleNames, function(sample) {
    # Skip empty sample names (e.g., empty lines in the input file)
    if (sample == "") {
      print("ignored_empty_lines")
    } else {
      # Path to the annotation file (gene model)
      anno.loc <- paste0(wd, "/1_intron_preprocessing/2_annotation_object/hg38_annotations_cds.rds")
      
      # If the annotation file is missing, print a warning
      if (!file.exists(anno.loc)) {
        print("annotation file is missing!")
      }
      
      # Load the annotation RDS file
      rdsSource <- readRDS(anno.loc)
      
      # Path to the RPKM expression data for the current sample
      exprs.loc <- paste0(
        wd, "/2_gene_preprocessing/3_gene_expression/results/", sample.name,
        "/gene_expression/", sample, "_rpkm.rds"
      )
      
      # Load the expression data (SummarizedExperiment object)
      exprGenes <- readRDS(exprs.loc)
      
      # Create a data.table for counts and RPKM values
      expr_count <- data.table()
      expr_count$count <- assays(exprGenes)$count
      expr_count$rpkm <- assays(exprGenes)$rpkm
      
      # Extract gene annotation information as a data.table
      expr_id_data <- as.data.table(rowData(exprGenes))
      
      # Combine gene annotation with expression data
      gene_exprs_data <- cbind(expr_id_data, expr_count)
      
      # Generate column names for count and RPKM columns
      count_nam <- paste(colnames(exprGenes), "_count", sep = "")
      rpkm_nam <- paste(colnames(exprGenes), "_rpkm", sep = "")
      
      # Assign new column names to the combined data
      names(gene_exprs_data) <- c("entrez.id", "symbol", count_nam, rpkm_nam)
      
      # Select only the entrez.id and count columns
      gene_exprs_count <- gene_exprs_data %>% select(entrez.id, count_nam)
      
      # Select only the entrez.id and RPKM columns
      gene_exprs_rpkm <- gene_exprs_data %>% select(entrez.id, rpkm_nam)
      
      # Path to the filtered intron annotation file
      rnintrons.file.loc <- paste0(
        wd, "/1_intron_preprocessing/3_filtering_gobj/rnhg38_filtered_introns_cds.rds"
      )
      
      # Load the filtered intron annotation data
      rn_introns <- readRDS(rnintrons.file.loc)
      
      # Convert the intron annotation to a data.frame and add a rowname column
      se_rn_introns <- as.data.frame(rn_introns)
      se_rn_introns$id <- rownames(se_rn_introns)
      
      # Select only the entrez.id and id columns, and order them
      se_rn_introns <- se_rn_introns %>% select(entrez.id, id)
      se_rn_introns <- se_rn_introns[, c(head(names(se_rn_introns), 1), tail(names(se_rn_introns), 1))]
      
      # Merge the intron annotation with the gene expression counts by entrez.id
      count <- merge(se_rn_introns, gene_exprs_count, by = "entrez.id", all = TRUE)
      count <- count[!is.na(count$id), ]  # Keep only rows with valid intron IDs
      rownames(count) <- count$id         # Set rownames to intron IDs
      count <- count %>% select(all_of(count_nam))  # Keep only the count columns
      
      # Merge the intron annotation with the gene expression RPKM by entrez.id
      rpkm <- merge(se_rn_introns, gene_exprs_rpkm, by = "entrez.id", all = TRUE)
      rpkm <- rpkm[!is.na(rpkm$id), ]    # Keep only rows with valid intron IDs
      rownames(rpkm) <- rpkm$id          # Set rownames to intron IDs
      rpkm <- rpkm %>% select(all_of(rpkm_nam))  # Keep only the RPKM columns
      
      # Define the output directory for results
      results.sample.dir <- paste0(gene_expression_data.dir, "/")
      
      # Save the processed count data as an RDS file
      saveRDS(count, paste0(results.sample.dir, count_nam, "_df", ".rds"))
      
      # Save the processed RPKM data as an RDS file
      saveRDS(rpkm, paste0(results.sample.dir, rpkm_nam, "_df", ".rds"))
      
      # Print a success message
      print("successfully retrived data")
    }
  })
}


geneexpr_se <- function(input.data.path, wd, sample.name) {
  # Function to combine gene expression count and RPKM data for a sample group
  # into SummarizedExperiment objects and save them for downstream analysis

  # Define the directory containing individual gene expression result files
  gene_expression_data.dir <- file.path(
    paste0(wd, "/2_gene_preprocessing/3_gene_expression/results/", sample.name, "/gene_expression_results")
  )
  
  # Define the directory to save combined SummarizedExperiment results
  gene_expression_results.dir <- file.path(
    paste0(wd, "/2_gene_preprocessing/3_gene_expression/results/", sample.name, "/gene_expression_se_results")
  )
  
  # Set output directory paths
  results.se.dir <- paste0(gene_expression_results.dir, "/")
  file_path <- paste0(gene_expression_data.dir, "/")
  
  # ----------------------------------------
  # Combine all *_count_df.rds files into a single data matrix
  # ----------------------------------------
  count <- do.call(
    cbind, 
    lapply(
      list.files(path = file_path, pattern = "*_count_df.rds", full.names = TRUE), 
      readRDS
    )
  )
  
  # Order the rows using mixedorder (from gtools) for consistency
  count <- count[mixedorder(rownames(count)), ]
  
  # Replace any NA values with "not_expressed"
  count[is.na(count)] <- "not_expressed"
  
  # Create a SummarizedExperiment object for counts
  se_count = SummarizedExperiment(assays = list(count = count))
  
  # Save the SummarizedExperiment object and the count matrix as RDS and CSV
  saveRDS(se_count, paste0(file_path, "se_gene_expr_count", ".RDS"))
  write.csv(count, paste0(file_path, "se_gene_expr_count", ".csv"))
  
  # ----------------------------------------
  # Combine all *_rpkm_df.rds files into a single data matrix
  # ----------------------------------------
  rpkm <- do.call(
    cbind, 
    lapply(
      list.files(path = file_path, pattern = "*_rpkm_df.rds", full.names = TRUE), 
      readRDS
    )
  )
  
  # Order the rows using mixedorder for consistency
  rpkm <- rpkm[mixedorder(rownames(rpkm)), ]
  
  # Replace any NA values with "not_expressed"
  rpkm[is.na(rpkm)] <- "not_expressed"
  
  # Create a SummarizedExperiment object for RPKM values
  se_rpkm = SummarizedExperiment(assays = list(rpkm = rpkm))
  
  # Save the SummarizedExperiment object and the RPKM matrix as RDS and CSV
  saveRDS(se_rpkm, paste0(results.se.dir, "se_gene_expr_rpkm", ".RDS"))
  write.csv(rpkm, paste0(results.se.dir, "se_gene_expr_rpkm", ".csv"))
  
  # Print a success message
  print("successfully created gene expression summarized experiment")
}
