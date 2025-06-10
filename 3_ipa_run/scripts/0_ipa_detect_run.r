

####################################################################################################################
# Set project directory
####################################################################################################################

project.dir <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"
setwd(project.dir)

####################################################################################################################
# Define the Data table
####################################################################################################################


# dt.name <- "data_table_test.txt"
# atlas_name <- "test"

# dt.name <- "data_table_test2.txt"
# atlas_name <- "test2"

dt.name <- "data_table_CoMMpass_1.txt"
atlas_name <- "CoMMpass"

datainfo.location <- file.path(project.dir, "input_data_tables", dt.name)

####################################################################################################################
# Get the genomic coverage and intron retention
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/1_ipa_detect.r")

ipa_detect(datainfo.location,project.dir, atlas_name)

retrieve_intronreten_data(datainfo.location,project.dir, atlas_name)

intronret_se(project.dir, atlas_name)

####################################################################################################################
# Get the introns with some coveage (5 reads for 100bp contiguous stretch)
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/2_filtering.r")
filtering(datainfo.location,project.dir, atlas_name)


####################################################################################################################
# Detect changepoints using PELT
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/3_pelt.r")
pelt(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# filter changepoints to detect IPAs
####################################################################################################################


source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/4_filter_cpts_de.R")
filter_changepoints(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Merge IPAs for all the chromosomes
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/5_merge.r")
merge_cpts(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Get the TPM clculation for new terminal exon
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/6_analyse_exon_structure.r")
filter_te_run(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Make IPA atlas per sample group
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/7_make_atlas.r")
make_atlas_run(datainfo.location,project.dir, atlas_name)


####################################################################################################################
# Get the IPA usage in all the samples
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/8_calculate_ipa_usage.R")
calc_ipa_usage(datainfo.location,project.dir, atlas_name)

####################################################################################################################
# Get SummarizedExperiment object for IPA usage accross the samples
####################################################################################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/3_ipa_run/scripts/9_create_ipa_usage_se.R")
ipa_usage_se(datainfo.location,project.dir, atlas_name)






