
####################################################################################################################
# Set project directory
####################################################################################################################

project.dir <- "/IPAseek_pipeline"
setwd(project.dir)

####################################################################################################################
# Define the Data table
####################################################################################################################

dt.name <- "data_table_test_uniq.txt"
atlas_name <- "test"

# dt.name <- "data_table_test2.txt"
# atlas_name <- "test2"


datainfo.location <- file.path(project.dir, "input_data_tables", dt.name)

####################################################################################################################
# Get the genomic coverage and intron retention
####################################################################################################################

source("/IPAseek_pipeline/3_ipa_run/scripts/1_ipa_detect.r")

## Step 1
ipa_detect(datainfo.location,project.dir, atlas_name)

## Step 2
retrieve_intronreten_data(datainfo.location,project.dir, atlas_name)

## Step 3
intronret_se(project.dir, atlas_name)

####################################################################################################################
# Get the introns with some coveage (5 reads for 100bp contiguous stretch)
####################################################################################################################

## Step 4
source("/IPAseek_pipeline/3_ipa_run/scripts/2_filtering.r")
filtering(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Detect changepoints using PELT
####################################################################################################################

## Step 5
source("/IPAseek_pipeline/3_ipa_run/scripts/3_pelt.r")
pelt(datainfo.location,project.dir, atlas_name)



####################################################################################################################
# filter changepoints to detect IPAs
####################################################################################################################

## Step 6
source("/IPAseek_pipeline/3_ipa_run/scripts/4_filter_cpts_de.R")
filter_changepoints(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Merge IPAs for all the chromosomes
####################################################################################################################

## Step 7
source("/IPAseek_pipeline/3_ipa_run/scripts/5_merge.r")
merge_cpts(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Get the TPM clculation for new terminal exon
####################################################################################################################

## Step 8
source("/IPAseek_pipeline/3_ipa_run/scripts/6_analyse_exon_structure.r")
filter_te_run(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Make IPA atlas per sample group
####################################################################################################################

## Step 9
source("/IPAseek_pipeline/3_ipa_run/scripts/7_make_atlas.r")
make_atlas_run(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Get the IPA usage in all the samples
####################################################################################################################

# ## Step 10
source("/IPAseek_pipeline/3_ipa_run/scripts/8_calculate_ipa_usage_combined.R")
calc_ipa_usage(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Get SummarizedExperiment object for IPA usage accross the samples
####################################################################################################################

## Step 11
source("/IPAseek_pipeline/3_ipa_run/scripts/9_create_ipa_usage_se.R")
ipa_usage_se(datainfo.location,project.dir, atlas_name)

# ####################################################################################################################
# # Calculate performace for test datasets
# ####################################################################################################################
# 
source("/IPAseek_pipeline/3_ipa_run/scripts/cal_precision.r")
cal_precision(datainfo.location,project.dir, atlas_name)
# 


# make all PROJECT_DIR="/IPAseek_pipeline" DATA_TABLE="/IPAseek_pipeline/input_data_tables/data_table_test2.txt" ATLAS_NAME="test2"

