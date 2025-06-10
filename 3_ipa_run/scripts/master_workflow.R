# master_workflow.R

### Job Tracking Functions ###
submit_jobs <- function(step_function = ipa_detect, datainfo.location = datainfo.location, project.dir = project.dir, atlas_name = atlas_name) {

  # Run the step function (which submits jobs internally)
  step_function(datainfo.location,project.dir, atlas_name)
  
  # Capture all job IDs for this step (modify functions to write to a file)
  job_ids <- readLines("current_jobs.log")
  return(job_ids)
}

wait_for_jobs <- function(job_ids) {
  while(TRUE) {
    # Check job status using sacct
    status <- system( paste("sacct -j", paste(job_ids, collapse=","), "--format=State -P -n"))
    
    # Exit loop when all jobs are completed
    if(all(status %in% c("COMPLETED", "FAILED", "CANCELLED"))) break
    
    # Wait before checking again
    Sys.sleep(60)
  }
}


### Main Execution ###
project.dir <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"
datainfo.location <- file.path(project.dir, "input_data_tables", "data_table_test.txt")
atlas_name <- "test"

# Step 1: IPA Detection & Intron retention Calculation
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/1_ipa_detect.r")
job_ids <- submit_jobs(ipa_detect, datainfo.location, project.dir, atlas_name)
wait_for_jobs(job_ids)
file.remove("current_jobs.log")

retrieve_intronreten_data(datainfo.location,project.dir, atlas_name)
intronret_se(project.dir, atlas_name)

# Step 2: Get the introns with some coveage (5 reads for 100bp contiguous stretch)
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/2_filtering.r")
job_ids <- submit_jobs(filtering, datainfo.location,project.dir, atlas_name)
wait_for_jobs(job_ids)

# Step 3: Detect changepoints using PELT
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/3_pelt.r")
job_ids <- submit_jobs(pelt,datainfo.location,project.dir, atlas_name)
wait_for_jobs(job_ids)

# Step 4: filter changepoints to detect IPAs
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/4_filter_cpts_de.R")
job_ids <- submit_jobs(filter_changepoints, datainfo.location,project.dir, atlas_name)
wait_for_jobs(job_ids)

# Step 5: Merge IPAs for all the chromosomes
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/5_merge.r")
merge_cpts(datainfo.location,project.dir, atlas_name)

# Step 6: Get the TPM clculation for new terminal exon
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/6_analyse_exon_structure.r")
job_ids <- submit_jobs(filter_te_run, datainfo.location,project.dir, atlas_name)
wait_for_jobs(job_ids)

# Step 7: Make IPA atlas per sample group
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/7_make_atlas.r")
job_ids <- submit_jobs(make_atlas_run, datainfo.location,project.dir, atlas_name)
wait_for_jobs(job_ids)

# Step 8: Get the IPA usage in all the samples
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/8_calculate_ipa_usage.R")
job_ids <- submit_jobs(calc_ipa_usage, datainfo.location,project.dir, atlas_name)
wait_for_jobs(job_ids)

# Step 9: Get SummarizedExperiment object for IPA usage accross the samples
source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/9_create_ipa_usage_se.R")
ipa_usage_se(datainfo.location,project.dir, atlas_name)

