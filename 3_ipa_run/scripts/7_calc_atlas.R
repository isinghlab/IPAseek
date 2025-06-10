
library(dplyr)          
library(data.table)     
library(GenomicRanges)  
library(splicejam)      

#################################################################################################################################  
##########################   Calculate IPA atlas accross multiple samples   #####################################################
#################################################################################################################################  

calc_atlas <- function(data.input){

   # Extract key parameters from the input data frame
   data.input_name <- data.input$atlas[1]         # Atlas batch name
   wd <- data.input$wd[1]                         # Working directory
   atlas_name <- data.input$atlas_name[1]         # Atlas name identifier

   # Source external R script with required functions for CDS overlap calculation
   source(paste0(wd,"/3_ipa_run/scripts/7_calc_cds_overlap.R"))

   # Aggregate sample names for each condition into a comma-separated string
   groups_df <- aggregate(NAME ~ CONDITION, data.input[c("CONDITION", "NAME")], function(x) paste(x, collapse = ","))

   # Extract sample names as a character vector
   sampleNames <- as.character(data.input$NAME)

   print(paste("processing for :", sampleNames)) # Track progress

   # Read and merge all exon expression CSVs for all samples in the atlas
   cpt_all <- do.call(
     rbind,
     lapply(
       list.files(
         path = paste0(wd,"/pelt/results/",atlas_name,"/",sampleNames,"/exon_exprs_results/"),
         full.names = TRUE,
         pattern = "*.csv"
       ),
       read.csv
     )
   )
   cpt_all <- unique(cpt_all) # Remove duplicate rows
   # cpt_all <- cpt_all[-1]   # (Optionally remove first row if needed)
   # Select relevant columns for downstream analysis
   cpt_all <- cpt_all[c("seqnames", "start", "end", "strand", "entrez.id", "id", "ipa_sel", "source")]

   # Map each terminal exon event to its experimental condition using sample name
   cpt_all$source_group <- data.input$CONDITION[match(cpt_all$source, data.input$NAME)]

   # Create a GRanges object for all terminal exons
   te_gr <- makeGRangesFromDataFrame(cpt_all, keep.extra.columns=TRUE)

   # Convert entrez.id to character and create a unique split column for each IPA event
   te_gr$entrez.id <- as.character(te_gr$entrez.id)
   te_gr$split_col <- paste(te_gr$id, te_gr$ipa_sel, sep=" ")

   # Split the GRanges object by unique IPA event (id + ipa_sel)
   te_gr_split <- split(te_gr, f = c(te_gr$split_col))

   # Optionally, you could limit to a subset for testing
   # te_gr_split <- te_gr_split[1:50]

   # Reduce IPA events within 100bp to a single representative region for the atlas
   print(paste("Making Atlas"))

   library(parallel) # Load parallel library for mclapply

   te_median <- mclapply(te_gr_split, function(gr){

      gr_orig <- gr                        # Store the original group
      gr <- reduce(gr, min.gapwidth=100L)  # Merge regions within 100bp

      gr <- annotateGRfromGR(gr, gr_orig)  # Annotate merged region with original metadata

      # For each merged region, set the terminal exon boundary to the median of the original ends
      gr_split <- split(gr)
      gr_median <- lapply(gr_split, function(gr_splitted){

        # For positive strand, use the median of the ends
        if(as.character(strand(gr_splitted)) == "+" ){
          gr_within <- subsetByOverlaps(gr_orig, gr_splitted)
          ends <- c(end(gr_within))
          ends_median <- median(ends)
          end(gr_splitted) <- ends_median
        } else { # For negative strand, use the median of the starts
          gr_within <- subsetByOverlaps(gr_orig, gr_splitted)
          ends <- c(start(gr_within))
          ends_median <- median(ends)
          start(gr_splitted) <- ends_median
        }
      })
      print(gr) # Print progress for each group
   })

   # Calculate retained CDS for each merged region using a custom function
   print("calculating retained cds")
   te_median <- mclapply(te_median, calc_cds_overlap)

   # Convert the list of GRanges to a single GRanges object
   te_median <- GRangesList(te_median)
   te_median <- unlist(te_median)

   # Count the number of unique sources (samples) supporting each atlas event
   te_median$no_of_sources <- sapply(strsplit(te_median$source,','), uniqueN)
   # Calculate the percentage of samples supporting each event
   te_median$percentage_of_sources <- (te_median$no_of_sources / length(sampleNames)) * 100

   # Calculate confidence for each atlas event using a custom function and group info
   te_median_confident <- calc_confidence(te_median, groups_df)

   # Prepare the output directory for the atlas if it doesn't exist
   res.dir <- file.path(wd, "pelt", "results", atlas_name)
   if(!dir.exists(res.dir)){
      dir.create(res.dir)
   }

   print(paste("Saving Atlas"))
   print(res.dir)
   print(atlas_name)
   print(data.input_name)

   # Save the full and confident atlas as RDS and CSV files for downstream use
   saveRDS(te_median, paste0(res.dir, "/",atlas_name,"_",data.input_name, "_ipa_atlas.RDS"))
   write.csv(as.data.frame(te_median), paste0(res.dir, "/",atlas_name,"_",data.input_name,"_ipa_atlas.csv"))

   saveRDS(te_median_confident, paste0(res.dir, "/",atlas_name,"_",data.input_name, "_ipa_atlas_conf.RDS"))
   write.csv(as.data.frame(te_median_confident), paste0(res.dir, "/",atlas_name,"_",data.input_name,"_ipa_atlas_conf.csv"))
}
