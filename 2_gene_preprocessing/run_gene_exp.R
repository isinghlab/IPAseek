# ── Parse command-line arguments ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

project.dir      <- parse_arg(args, "--project_dir",
                               default = Sys.getenv("IPASEEK_PROJECT_DIR",
                                         unset = normalizePath(".")))
dt.name          <- parse_arg(args, "--data_table",      default = "input_data_tables/data_table_test.txt")
sample.name      <- parse_arg(args, "--atlas_name",      default = "ipaseek_run")
star.index.path  <- parse_arg(args, "--star_index",
                               default = Sys.getenv("IPASEEK_STAR_INDEX", unset = ""))
annotation.dir.arg <- parse_arg(args, "--annotation_dir",
                               default = Sys.getenv("IPASEEK_ANNOTATION_DIR", unset = ""))
# ─────────────────────────────────────────────────────────────────────────────

# Set env vars so sourced scripts pick them up
if (nchar(project.dir) > 0)        Sys.setenv(IPASEEK_PROJECT_DIR    = project.dir)
if (nchar(star.index.path) > 0)    Sys.setenv(IPASEEK_STAR_INDEX     = star.index.path)
if (nchar(annotation.dir.arg) > 0) Sys.setenv(IPASEEK_ANNOTATION_DIR = annotation.dir.arg)

#########################################################
## Get data_tables
#########################################################

# If data_table is a relative path, make it absolute relative to project.dir
if (!file.exists(dt.name)) {
  dt.name <- file.path(project.dir, dt.name)
}
datainfo.location <- dt.name

#########################################################
## 1 -- Align fastqs to the reference genome
#########################################################

source(file.path(project.dir, "2_gene_preprocessing", "1_rnaseq_pipeline", "rnaseq_pipeline.R"))

pushSTAR(datainfo.location)

#########################################################
## 2 -- Create unique bams
#########################################################

source(file.path(project.dir, "2_gene_preprocessing", "2_bams", "scripts", "filt_bam_uniq.R"))

setwd(project.dir)

filt_bam(datainfo.location,project.dir)

# ── Auto-generate the _uniq data table after BAM filtering ───────────────────
generate_uniq_table <- function(dt.path, project.dir) {
  dt <- read.delim(dt.path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  # Update FILE_PATH to point to the uniq_bams subdirectory
  dt$FILE_PATH <- file.path(dt$FILE_PATH, "uniq_bams")
  # Compute LIB_SIZE from uniquely mapped BAM if samtools is available
  dt$LIB_SIZE <- sapply(seq_len(nrow(dt)), function(i) {
    bam <- file.path(dt$FILE_PATH[i], paste0(dt$NAME[i], "_uniq.bam"))
    if (file.exists(bam) && nchar(Sys.which("samtools")) > 0) {
      as.integer(system(paste("samtools view -c -F 4", bam), intern=TRUE))
    } else {
      NA_integer_
    }
  })
  # Write the uniq table next to the original
  uniq.path <- sub("\\.txt$", "_uniq.txt", dt.path)
  write.table(dt, uniq.path, sep="\t", quote=FALSE, row.names=FALSE)
  message(sprintf("Auto-generated uniq data table: %s", uniq.path))
  return(uniq.path)
}

datainfo.location_uniq <- generate_uniq_table(datainfo.location, project.dir)

#########################################################
## 3 -- Get RNA-seq counts
#########################################################

source(file.path(project.dir, "2_gene_preprocessing", "1_rnaseq_pipeline", "rnaseq_pipeline.R"))

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


source(file.path(project.dir, "2_gene_preprocessing", "1_rnaseq_pipeline", "gene_expr.R"))

start_time <- Sys.time()
setwd(project.dir)
count.files.dir <- file.path(paste0(project.dir,"/2_gene_preprocessing/3_gene_expression/results/",sample.name,"/rnaseq_counts"))

df.sample.info <- read.delim(datainfo.location_uniq, sep = "\t", stringsAsFactors = F)


se.gene <- createSE(project.dir, count.files.dir, data.table.path = datainfo.location_uniq, design.colmns = c("CONDITION"))
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


source(file.path(project.dir, "2_gene_preprocessing", "1_rnaseq_pipeline", "expressed_genes_data_retrieved.R"))

retrieve_geneexpr_data(datainfo.location_uniq,project.dir, sample.name)
geneexpr_se(project.dir,project.dir, sample.name)

