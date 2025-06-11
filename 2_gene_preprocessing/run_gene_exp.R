
#########################################################
## Define project directory
#########################################################

# project.dir <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"
project.dir <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"

#########################################################
## Get data_tables
#########################################################

dt.name <- "data_table_test.txt"
sample.name <- "test"

datainfo.location <- file.path(project.dir, "data_tables", dt.name)

#########################################################
## 1 -- Align fastqs to the reference genome
#########################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/rnaseq_pipeline.R")

pushSTAR(datainfo.location)

#########################################################
## 2 -- Create unique bams
#########################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/2_gene_preprocessing/2_bams/scripts/filt_bam_uniq.R")

setwd(project.dir)

filt_bam(datainfo.location,project.dir)

#########################################################
## Define project directory
#########################################################

project.dir <- "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline"

#########################################################
## Get uniq data_tables
#########################################################

# dt.name_uniq <- "data_table_test.txt"
# sample.name <- "test"

dt.name_uniq <- "data_table_CoMMpass_1.txt"
sample.name <- "CoMMpass"

datainfo.location_uniq <- file.path(project.dir, "input_data_tables", dt.name_uniq)

#########################################################
## 3 -- Get RNA-seq counts
#########################################################

source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/rnaseq_pipeline.R")

print(datainfo.location_uniq)


sample.dir <- file.path(paste0(project.dir,"/2_gene_preprocessing/3_gene_expression/results/",sample.name))

if(!dir.exists(sample.dir)){
 		dir.create(sample.dir)
 	}

count.obj.dir <- file.path(paste0(project.dir,"/2_gene_preprocessing/3_gene_expression/results/",sample.name,"/rnaseq_counts"))

if(!dir.exists(count.obj.dir)){
 		dir.create(count.obj.dir)
 	}


runCDSCounts(project.dir, datainfo.location_uniq, count.obj.dir)


#########################################################
## 4 -- Get SE objects for raw counts and RPKM
#########################################################


source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/gene_expr.R")    

start_time <- Sys.time()
setwd(project.dir)
count.files.dir <- file.path(paste0(project.dir,"/2_gene_preprocessing/3_gene_expression/results/",sample.name,"/rnaseq_counts"))

df.sample.info <- read.delim(datainfo.location_uniq, sep = "\t", stringsAsFactors = F)


se.gene <- createSE(project.dir, count.files.dir, data.table.path = datainfo.location_uniq, design.colmns = c("CONDITION", "RACE", "GENDER", "ETHNICITY", "VITAL_STATUS", "AGE_AT_INDEX", "ISS_STAGE", "TUMOR_DESCRIPTOR"))
se.gene <- getGeneSE(se.gene)

samples.of.interest <- df.sample.info$NAME 
exprs.files.dir <- file.path(paste0(project.dir, "/2_gene_preprocessing/3_gene_expression/results/",sample.name,"/gene_expression"))

if(is.null(exprs.files.dir)){
		exprs.files.dir <- file.path(paste0(project.dir,"/2_gene_preprocessing/3_gene_expression/results/",sample.name,"/gene_expression"))
	}
	
	if(!dir.exists(exprs.files.dir)){
 		dir.create(exprs.files.dir)
 	}


for(i in 1:length(samples.of.interest)){

se.gene.sample <- se.gene[ ,samples.of.interest[i]]
saveRDS(se.gene.sample,file=file.path(exprs.files.dir,paste0(samples.of.interest[i],"_all.rds")))
se.gene.sample.expr  <- getExpressedGenes(se.gene.sample, rpkm.cutoff = 0.5)
saveRDS(se.gene.sample.expr ,file=file.path(exprs.files.dir,paste0(samples.of.interest[i],"_rpkm.rds")))

}

end_time <- Sys.time()
exec.time<-end_time - start_time
print(paste0("execution time:",exec.time))

#########################################################
## 5 -- Retreive exp values from SE objects
#########################################################


source("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/expressed_genes_data_retrieved.R")

retrieve_geneexpr_data(datainfo.location_uniq,project.dir, sample.name)
geneexpr_se(project.dir,project.dir, sample.name)

