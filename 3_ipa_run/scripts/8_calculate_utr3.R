
library(GenomicAlignments)   
library(GenomicFeatures)     
library(dplyr)               
library(tidyr)               
library(data.table)          
library(gtools)              
library(ggplot2)            
library(scales)             

#################################################################################################################################  
#################   Function to get last CDS start/end coordinates for gene with IPA event ######################################
#################################################################################################################################  

calc_cds <- function(ipa, wd){

   print("calc_cds")             

   # ipa: a GRanges object representing a single IPA event
   # wd: working directory path

   # Extract the entrez gene ID from the IPA event
   ipa_entrezid <- ipa$entrez.id
   print(ipa_entrezid)          # Print the gene ID for tracking

   # Load the reference annotation for the genome (must contain CDS/exon info)
   hg38 <- readRDS(paste0(wd, "/1_intron_preprocessing/1_flatten_genome/hg38_annotated_numbered_cds.rds"))

   # Subset the annotation to only this gene
   hg38_ipa <- hg38[hg38$entrez.id %in% ipa_entrezid,]

   # Check if there is CDS annotation for this gene
   if("cds" %in% hg38_ipa$exon.anno){

      # Extract only CDS features for this gene
      hg38_ipa_cds <- hg38_ipa[hg38_ipa$exon.anno %in% c("cds")]
      # Take the last CDS entry (typically the most 3' CDS, relevant for IPA)
      hg38_ipa_cds <- tail(hg38_ipa_cds, 1)

      # Calculate the width of the IPA and CDS regions (not used further here)
      width_ipa <- width(ipa)
      width_cds <- width(hg38_ipa_cds)

      # Annotate the IPA event with the CDS start and end coordinates
      ipa$cds_start <- start(hg38_ipa_cds)
      ipa$cds_end <- end(hg38_ipa_cds)

      # Optionally, print widths for debugging
      # print(paste("width ipa : ", ipa_width, "| width utr3 :", utr3_width))
   } else {
      # If no CDS annotation is present, set CDS start/end to "none"
      ipa$cds_start <- "none"
      ipa$cds_end <- "none"
   }

   # Return the annotated IPA object
   ipa
}

