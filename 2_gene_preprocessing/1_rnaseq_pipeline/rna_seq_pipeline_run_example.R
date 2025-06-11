source("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/rnaseq_pipeline.R")

##############
project.dir <- "/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/"

########

# #dt.name <- "data_table_GSE96800.txt"
# #dt.name <- "data_table_GSE213909.txt"
# dt.name <- "data_table_GSE125534.txt"
# #dt.name <- "data_table_GSE140528.txt"
# #dt.name <- "data_table_GSE138913.txt"


#dt.name <- "data_table_GSE173790.txt"

#dt.name <- "data_table_GSE79044.txt"

#dt.name <- "encode_rnaseq.txt"


#dt.name <- "data_table_encode8.txt"


datainfo.location <- file.path(project.dir, "data_tables", dt.name)
pushSTAR(datainfo.location)
#count.obj.dir <- "/storage/cylin/grail/projects/tamu/09182019/rnaseq_counts"
#runCDSCounts(project.dir, datainfo.location, count.obj.dir)
