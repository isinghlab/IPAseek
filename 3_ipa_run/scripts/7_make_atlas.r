
library(dplyr)          
library(data.table)     
library(GenomicRanges)   
library(splicejam)       

# Function to set up and submit a SLURM job for running the make_atlas function
make_atlas_run <- function(input.data.path, wd, atlas_name){

  # Define paths for SLURM submission scripts and log files
  slurm.files.dir <- file.path(wd, "pelt", "slurm_submission_make_atlas", atlas_name)
  logs.files.dir  <- file.path(wd, "pelt", "logs_make_atlas", atlas_name)

  # Create SLURM and logs directories if they do not already exist
  if(!dir.exists(slurm.files.dir)){
    dir.create(slurm.files.dir, recursive = TRUE)
  }
  if(!dir.exists(logs.files.dir)){
    dir.create(logs.files.dir, recursive = TRUE)
  }

  # Read the input metadata file (tab-delimited), which contains sample names and BAM file locations
  data.input <- read.delim(input.data.path, sep="\t", header=TRUE)

  # Prepare the paths object with necessary metadata for the atlas run
  paths <- data.input
  paths$atlas_name <- atlas_name
  paths$wd <- wd

  # Save the make_atlas function environment for later use in the SLURM job
  make_atlas.loc <- paste0(slurm.files.dir, "/make_atlas.Rdata")
  save("make_atlas", file=make_atlas.loc)

  # Create a sample-specific directory for this atlas if it does not exist
  if(!dir.exists(paste0(slurm.files.dir, "/", atlas_name))){
    dir.create(paste0(slurm.files.dir, "/", atlas_name))
  }

  # Define paths for the R script and the sample directory
  sample.dir  <- file.path(paste0(slurm.files.dir, "/", atlas_name, "/"))
  script.name <- file.path(paste0(sample.dir, '/', "_make_atlas.R"))

  # Save the paths object for use in the SLURM job
  save(paths, file=paste0(sample.dir, atlas_name, "_", "paths.Rdata"))

  # Begin writing the R script that will be executed by SLURM
  sink(file=script.name)
  # Write library loading commands to the script
  cat("
      \nlibrary(dplyr)
      \nlibrary(data.table)
      \nlibrary(GenomicRanges)
      \nlibrary(splicejam)")
  # Load the make_atlas function environment
  cat(paste0("\nload ('", make_atlas.loc, "')"))
  # Load the previously saved paths object
  cat(paste0("\nload ('", sample.dir, atlas_name, "_", "paths.Rdata')"))
  # Call the make_atlas function with the loaded paths object
  cat(paste0("\nmake_atlas(paths)"))
  sink() # Close the script file

  # Create the SLURM submission bash script
  bash.file.location <- file.path(paste0(sample.dir, 'submit.sh'))
  slurm.jobname <- sprintf("#SBATCH --job-name=%s_make_atlas", atlas_name)      
  slurm.time    <- sprintf("#SBATCH --time=36:00:00")                           
  slurm.mem     <- sprintf("#SBATCH --mem=48G")                                  
  slurm.tasks   <- sprintf("#SBATCH --ntasks=8")                                
  slurm.output  <- paste0("#SBATCH --output=", logs.files.dir, "/", atlas_name, "_make_atlas") 
  sbatch.line   <- sprintf("Rscript --vanilla %s", script.name)              

  # Write all SLURM directives and the Rscript command to the bash file
  file.conn <- file(bash.file.location)
  writeLines(
    c(
      "#!/bin/bash",
      slurm.jobname,
      slurm.time,
      slurm.mem,
      slurm.tasks,
      slurm.output,
      sbatch.line
    ),
    file.conn
  )
  close(file.conn)

  # Submit the SLURM job using the sbatch command
  system(paste0("sbatch ", bash.file.location))

  # Print a message indicating that the job was submitted
  print(paste0("Running ... ", "job for ", atlas_name, " submitted."))

}

#################################################################################################################################  
##########################   Function to call the creation of a IPA atlas   #####################################################
#################################################################################################################################  
 
make_atlas <- function(paths) {

    # Extract key parameters from the input paths object
    atlas_name <- paths$atlas_name[[1]]  # Name identifier for the atlas
    wd <- paths$wd[[1]]                  # Working directory path
    data.input <- paths                  # Assign input data to a working variable

    # Split input data by CONDITION column to process groups separately
    data.input.split <- split(data.input, data.input$CONDITION)
    # Add a "full" category containing the unsplit/combined dataset
    data.input.split$full <- data.input

    # Source external script containing the core atlas calculation logic
    # Note: This script must define the calc_atlas() function
    source(paste0(wd, "/3_ipa_run/scripts/7_calc_atlas.R"))

    # Assign atlas batch names to each split group
    for (i in 1:length(data.input.split)) {
        data.input.split[[i]]$atlas <- names(data.input.split)[i]
    }

    # Reassign to the "full" combined dataset for final processing
    data.input.split <- data.input.split$full

    # Execute the core atlas calculation function with prepared data
    calc_atlas(data.input.split)
}


