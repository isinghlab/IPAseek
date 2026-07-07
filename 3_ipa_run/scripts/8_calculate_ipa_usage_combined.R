
library(GenomicAlignments)
library(dplyr)
library(tidyr)
library(data.table)
library(gtools)
library(ggplot2)

      
calc_ipa_usage <- function(input.data.path, wd, atlas_name, mode = "slurm"){
# 
# 
# input.data.path <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/input_data_tables/data_table_combined_uniq.txt"
# wd <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"
# atlas_name <- "combined"
slurm.files.dir <- file.path(wd,"pelt", "slurm_submission_ipa_usage")
slurm.files.dir <- file.path(wd,"pelt", "slurm_submission_ipa_usage")
logs.files.dir <- file.path(wd,"pelt", "logs_ipa_usage")

if(!dir.exists(slurm.files.dir )){
if(!dir.exists(slurm.files.dir )){
        dir.create(slurm.files.dir )
   }
if(!dir.exists(logs.files.dir )){
if(!dir.exists(logs.files.dir )){
        dir.create(logs.files.dir )
   }


data.input <- read.delim(input.data.path, sep="\t", header=T)
data.input <- read.delim(input.data.path, sep="\t", header=T)
sampleNames <- paste0( as.character(data.input$ATLAS_NAME),"," ,as.character(data.input$NAME))
# sampleNames <- tail(sampleNames,15)

  #Run for each sample 
  sapply(sampleNames, function(sample){
    


    # sample = "encode,ENCFF831UDC"

      atlas_name1 = strsplit(sample, split = ",")[[1]][1]
      sample = strsplit(sample, split = ",")[[1]][2]
      
      res_file <- paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_", sample, "_ipa_usage_atlas.csv")
      if(!file.exists(res_file)){

      ## Get the data
      filePath <- data.input[data.input$NAME == sample, ]$FILE_PATH
      reads.dir <- file.path(filePath, paste0(sample, "_reads.rds"))
      
      ipa_atlas_path <-  paste0(wd,"/pelt/results/",atlas_name,"/", atlas_name,"_full_ipa_atlas_conf.RDS")
      
      #cpts_path<-paste0(wd,"/pelt/results/",sample,"/filter_cpts_results/")
      library_size <- read.table(paste0(wd ,"/pelt/results/",atlas_name1,"/", sample,"_library_size"))
      library_size <-  library_size$V1


      paths <- data.frame(
        reads.dir = reads.dir,
        sample = sample,
        condition = data.input[data.input$NAME == sample, "CONDITION"],
        ipa_atlas = ipa_atlas_path,
        wd = wd,
        ipa_atlas_name = atlas_name,
        ipa_atlas_name = atlas_name1,
        library_size = library_size
      )

      if (mode == "local") {
        ipa_usage(paths)
        print(paste0("Running ... IPA usage for ", sample, " (local mode)."))
      } else {
      ipa_usage.loc<-paste0(slurm.files.dir,"/ipa_usage.Rdata")
      save("ipa_usage", file=ipa_usage.loc)
          
      if(!dir.exists(paste0(slurm.files.dir,"/", sample))){
                dir.create(paste0(slurm.files.dir,"/", sample))
           }
      if(!dir.exists(logs.files.dir)){
                dir.create(logs.files.dir)
          }


      sample.dir <- file.path(paste0(slurm.files.dir,"/", sample, "/"))
      script.name <- file.path(paste0(sample.dir,'/', "_ipa_usage.R"))
      save(paths, file=paste0(sample.dir,sample,"_","paths.Rdata"))
      sink(file=script.name)

      cat(paste0("\nload (\'",ipa_usage.loc,"')"))
      #load the previously saved paths object
      cat(paste0("\nload (\'",sample.dir,sample,"_", "paths.Rdata\')"))
      #call the function
      cat(paste0("\nipa_usage(paths)"))
      #close
      sink()

      #create the bash script
            bash.file.location <- file.path(paste0(sample.dir,'submit.sh'))
            slurm.jobname <- sprintf("#SBATCH --job-name=%s_ipa_usage", sample)
            slurm.time <- sprintf("#SBATCH --time=15:00:00")
            slurm.mem <- sprintf("#SBATCH --mem=16G")
            slurm.tasks <- sprintf("#SBATCH --ntasks=1")
            slurm.output <- paste0("#SBATCH --output=",logs.files.dir,"/",sample,"_ipa_usage")
            sbatch.line<- sprintf("Rscript --vanilla %s", script.name)
            file.conn <- file(bash.file.location)
            writeLines(c("#!/bin/bash", slurm.jobname, slurm.time, slurm.mem, slurm.tasks, slurm.output, sbatch.line), file.conn)
            close(file.conn)
            system(paste0("sbatch ", bash.file.location))

         
       print(paste0("Running ... ", "job for ", sample, " submitted."))
       } # end slurm branch

}
})




}



ipa_usage <- function(paths){

    # Extract relevant parameters and file paths from the input 'paths' object
    reads_path      <- paths$reads.dir          # Path to sample's aligned reads (RDS)
    sample          <- paths$sample             # Sample name
    ipa_atlas_path  <- paths$ipa_atlas          # Path to confident IPA atlas (RDS)
    atlas_name      <- paths$ipa_atlas_name     # Atlas name identifier
    atlas_name1     <- paths$ipa_atlas_name1
    library_size    <- paths$library_size       # Library size (for RPKM normalization)
    wd              <- paths$wd                 # Working directory
    
    # Source external script with required functions for UTR3 and CDS calculations
    source(paste0(wd,"/3_ipa_run/scripts/8_calculate_utr3.R"))
    
    print(paste("calculating IPA usage : ", sample)) # Progress message
    
    # Load the confident IPA atlas for this sample
    print("Getting IPA_atlas")
    ipa_atlas <- readRDS(ipa_atlas_path)
    
    # Split the atlas into a list of GRanges objects (one per event)
    ipa_atlas <- split(ipa_atlas, as.factor(ipa_atlas))
    
    # ipa_atlas <- head(ipa_atlas, 100)
    
    # For each IPA event, calculate CDS information using the external function
    library(parallel)
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
    
    print(paste("RPKM calculation"))
    
    # Calculate RPKM for IPA and CDS regions
    ipa_atlas_df <- as.data.frame(ipa_atlas) %>%
      mutate(
        ipa_rpkm = (ipa_reads * 1e9) /(ipa_width * as.numeric(library_size)) ,
        cds_rpkm = (cds_reads * 1e9) /(cds_width * as.numeric(library_size)) 
      )
    
    
    print(paste("IPA usage"))
    # Calculate IPA usage as the fraction of IPA RPKM over total (IPA + CDS) RPKM
    ipa_atlas_df$cds_ipa_usage <- (ipa_atlas_df$ipa_rpkm / (ipa_atlas_df$ipa_rpkm + ipa_atlas_df$cds_rpkm))
    
    print(paste("Saving for : ", sample))
    # Save the final annotated data frame as a CSV file for downstream analysis
    write.csv(
      ipa_atlas_df,
      paste0(wd, "/pelt/results/", atlas_name, "/", atlas_name, "_", sample, "_ipa_usage_atlas.csv")
    )
    
    print("finished") # Indicate completion for this sample


}






