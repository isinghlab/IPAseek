
library(dplyr)            
library(data.table)       
library(GenomicRanges)    
library(splicejam)        

#################################################################################################################################  
########################### Function to calculate the fraction of coding sequence (CDS) retained ################################
########################### by a terminal exon event, relative to the annotated CDS for the gene  ###############################
#################################################################################################################################  


calc_cds_overlap <- function(terminal_exon_gr){

    # terminal_exon_gr: a GRanges object representing terminal exon(s) for a gene

    # Extract the entrez gene ID from the terminal exon
    te_entrezid <- terminal_exon_gr$entrez.id[1]

    # Load pre-annotated CDS and exon data for hg38 genome (must be available at this path)
    hg38 <- readRDS("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/1_intron_preprocessing/1_flatten_genome/hg38_annotated_numbered_cds.rds")
    
    # Extract CDS features for this gene
    hg38_cds <- hg38[hg38$exon.anno == "cds" & hg38$entrez.id %in% te_entrezid,]

    # Extract all exons for this gene
    hg38_te <- hg38[hg38$entrez.id %in% te_entrezid,]

    # Optionally: Extract 3' UTR exons for this gene (not used)
    # hg38_te_utr3 <- hg38_te[hg38_te$exon.anno2 %in% c("utr3")]
    # print(length(hg38_te_utr3))

    # Limit to CDS, intron, UTR3, and UTR5 features for the gene
    hg38_te <- hg38_te[hg38_te$exon.anno %in% c("cds", "intron", "utr3", "utr5")]
    # Merge overlapping/adjacent exons with a large gap tolerance (100kb)
    hg38_te <- reduce(hg38_te, min.gapwidth=100000L)
    # print(hg38_te)

    # For positive strand, set start of transcript to gene start; for negative, set end to gene end
    if(as.character(strand(terminal_exon_gr)[1]) == "+"){
      # For plus strand, set start position to gene start
      ipa_transcript <- terminal_exon_gr
      start(ipa_transcript) <- start(hg38_te)
    } else{
      # For minus strand, set end position to gene end
      ipa_transcript <- terminal_exon_gr
      end(ipa_transcript) <- end(hg38_te)
    }

    # Calculate the percent of CDS retained in the IPA transcript
    # (width of intersection between transcript and CDS, divided by total CDS width)
    terminal_exon_gr$perc_cds_retained <- (sum(width(intersect(ipa_transcript, hg38_cds))) / sum(width(hg38_cds))) * 100

    # Return the annotated GRanges object
    terminal_exon_gr
}

# Function to assign confidence to each IPA event based on group/sample support
calc_confidence <- function(ipa_atlas, groups_df){

  # ipa_atlas: GRanges object with IPA events
  # groups_df: data frame mapping CONDITION to comma-separated sample names

  # Get all unique group/condition names
  groups <- unique(groups_df$CONDITION) 

  # Split the IPA atlas into a list of GRanges (one per event)
  ipa_atlas_split <- split(ipa_atlas)

  # For each IPA event, count support in each group and assign confidence
  ipa_atlas_count <- lapply(ipa_atlas_split, function(gr){

    # gr: a GRanges object for a single IPA event

    mcols(gr)$confidence <- FALSE  # Initialize confidence as FALSE

    # print(gr)

    # For each group, count how many unique samples support this IPA event
    for(group in groups) {

      # print(group)

      # Get all sample names for this group
      group_cols <- unlist(strsplit(groups_df[groups_df$CONDITION == group, "NAME"], ","))
      # Get all sample names supporting this event
      source_cols <- unlist(strsplit(gr$source, ","))

      # Count how many group samples are present in the source/support list
      count_present <- length(intersect(group_cols, source_cols))

      # Store the count as a metadata column for this group
      mcols(gr)[[paste0(group, "_count")]] <- count_present

      # If at least two samples from this group support the event, or if already confident, set TRUE
      if(mcols(gr)[[paste0(group, "_count")]] >= 2 | mcols(gr)$confidence == TRUE ){
        mcols(gr)$confidence <- TRUE
      } else {
        mcols(gr)$confidence <- FALSE
      }
    }

    print(gr) # Print the GRanges object (for debugging/progress)

  })

  # Combine the list of annotated GRanges back into a single GRanges object
  ipa_atlas_count <- GRangesList(ipa_atlas_count)
  ipa_atlas_count <- unlist(ipa_atlas_count) 

  # Create a vector of group count column names (not used here, but may be useful downstream)
  filt_cols <- paste0(groups, "_count")

  # Return the annotated IPA atlas with confidence columns
  return(ipa_atlas_count)
}
