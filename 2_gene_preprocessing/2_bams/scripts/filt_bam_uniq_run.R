source("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/gene_preprocessing/bams/scripts/filt_bam_uniq.R")

project.dir <- "/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline"

setwd(project.dir)

#dt.name <- "data_table_test.txt"
#dt.name <- "data_table_hg38.txt"

#dt.name <- "data_table_GSE173790.txt"

#dt.name <- "data_table_GSE66117.txt"
#dt.name <- "data_table_GSE96800.txt"
#dt.name <- "data_table_GSE213909.txt"
dt.name <- "data_table_GSE125534.txt"
#dt.name <- "data_table_GSE140528.txt"
#dt.name <- "data_table_GSE138913.txt"

#dt.name <- "encode_rnaseq.txt"


#dt.name <- "data_table_encode.txt"
dt.location <- file.path(project.dir, "data_tables", dt.name)


filt_bam(dt.location,project.dir)

