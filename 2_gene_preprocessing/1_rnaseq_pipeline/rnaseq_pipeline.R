options(warn = -1)
library(data.table)
library(GenomicFeatures)
library(GenomicAlignments)
library(parallel)

# ######################################
# STAR
# Indexes on the shared data 
star.indexes <- c(hg19 = '/scratch/user/richa.rashmi.1202/genomes/homosapiens/ucsc/hg19/sequence/starindex',
	hg38 = '/scratch/user/richa.rashmi.1202/genomes/homosapiens/ucsc/hg38/sequence/starindex')

pushSTAR <- function(dt.location, ncores = 24, additional.args = NULL){

	# dt.location = "/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/data_tables/data_table_CLL11.txt"
	# ncores = 24
	# additional.args = NULL

	#dt.location = datainfo.location
	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(df.sample.info$NAME)
	sapply(sample.names, function(sample.name){

		# sample.name = "RNA-seq_CLL11"

		print(sample.name)
		
		##Get the genome
		genome.name <- as.character(subset(df.sample.info, NAME == sample.name)$GENOME)

		##location of the fastq file
		fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
		fastq.files <- unlist(strsplit(split = "::", fastq.file.location))
		##Find the endedness of the fastqs
		if(length(fastq.files) == 2){
			endedness <- "pe"
		} else if (length(fastq.files) == 1){
			endedness <- "se"
		}

		##Check if the files exists
		if(!all(file.exists(fastq.files)))
		stop("Error: Fastq files missing")

		##location of the bam file
		file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)
		bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == file.id)$FILE_PATH)
  		all.bam.files <- list.files(path = bam.file.dir, pattern = "*.bam$")
		bam.file.logical <- grepl(pattern = file.id, x = all.bam.files, fixed = T)
		##Check if the bam file already exists
		if(all(!bam.file.logical)){
			aligned <- FALSE
		} else {
			bam.file.location <- file.path(bam.file.dir, all.bam.files[bam.file.logical])
			if(file.exists(bam.file.location)){
				print(sprintf("BAM already exists for %s", sample.name))
			}
			aligned <- TRUE
		}
		
		if(!aligned){
			##Star cmd	
			star.cmd <- sprintf ("STAR --genomeLoad LoadAndRemove --genomeDir %s", star.indexes[genome.name])
	
 			#Finds all the Paired end files
    		star.cmd <- sprintf ("%s --readFilesIn %s", star.cmd, fastq.files[1])
			if (endedness == "pe")
			    star.cmd <- sprintf ("%s %s", star.cmd, fastq.files[2])
			
			#Unzipping
			star.cmd <- sprintf ("%s --readFilesCommand zcat", star.cmd)
			# Parallel processing
			star.cmd <- sprintf ("%s --runThreadN %d", star.cmd, ncores)
		
			# Additional args
			star.cmd <- sprintf ("%s --alignIntronMin 70 --alignIntronMax 100000 --outFilterMultimapScoreRange 0", star.cmd)
			star.cmd <- sprintf ("%s --outFilterMultimapNmax 10 --outFilterMismatchNmax 2 --outSAMattributes All", star.cmd)
			star.cmd <- sprintf ("%s --outSAMattrIHstart 0", star.cmd)
			star.cmd <- sprintf ("%s --outSJfilterOverhangMin 8 8 8 8 --outSJfilterCountUniqueMin 1 1 1 1 --outReadsUnmapped Fastx", star.cmd)
			star.cmd <- sprintf ("%s --outSAMstrandField intronMotif", star.cmd)
			star.cmd <- sprintf ("%s --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate", star.cmd)
			star.cmd <- sprintf ("%s --outWigType wiggle --outWigStrand Unstranded --outWigNorm None", star.cmd)
			star.cmd <- sprintf ("%s --limitBAMsortRAM 53687091200", star.cmd)
		
			if (!is.null (additional.args)) 
			    star.cmd <- sprintf ("%s %s", star.cmd, additional.args)
	
			output.file <- file.id
        	star.cmd <- sprintf ("%s > %s.bam", star.cmd, output.file)
		
    		##Index command
			index.cmd <- sprintf("samtools index %s.bam", output.file)
		
			##Get the wiggle and bedGraph files
			star.signal.cmd <- sprintf("STAR --runMode inputAlignmentsFromBAM --inputBAMfile %s.bam", output.file)
			star.wig.cmd <- sprintf("%s --outWigType wiggle --outWigStrand Unstranded", star.signal.cmd)
			star.bg.cmd <- sprintf("%s --outWigType bedGraph --outWigStrand Unstranded", star.signal.cmd)
		
			##Rename the wiggle files
			mv.wig.cmd1 <- sprintf("mv Signal.Unique.str1.out.wig %s.unique.wig", output.file)
			mv.wig.cmd2 <- sprintf("mv Signal.UniqueMultiple.str1.out.wig %s.unique.multiple.wig", output.file)
			mv.bg.cmd1 <- sprintf("mv Signal.Unique.str1.out.bg %s.unique.bg", output.file)
			zip.bg.cmd1 <- sprintf("gzip %s.unique.bg", output.file)
			mv.bg.cmd2 <- sprintf("mv Signal.UniqueMultiple.str1.out.bg %s.unique.multiple.bg", output.file)
			zip.bg.cmd2 <- sprintf("gzip %s.unique.multiple.bg", output.file)
			mv.bam.cmd <- sprintf("mv %s.bam %s", output.file, bam.file.dir)
			mv.bai.cmd <- sprintf("mv %s.bam.bai %s", output.file, bam.file.dir)
	
			##creat and change directory
			output.dir <- file.path(bam.file.dir, sample.name)
			dir.create(output.dir)
			cd.cmd <- sprintf("cd %s", output.dir)
		
			bash.file.location <- file.path(output.dir, paste0(output.file, ".sh"))
			file.conn <- file(bash.file.location)
    		writeLines(c('#!/bin/bash', cd.cmd, star.cmd, index.cmd, 
    			star.bg.cmd, mv.bg.cmd1, zip.bg.cmd1, mv.bg.cmd2, 
    			zip.bg.cmd2, mv.bam.cmd, mv.bai.cmd), file.conn)
			close(file.conn)
	
			system(sprintf("chmod u+x %s", bash.file.location))
			system(sprintf("sbatch --nodes=1 --ntasks-per-node=8 --mem=56G --time=08:00:00 %s", bash.file.location))			
		} else {
			return(NULL)
		}
	})
}

####################################################################################################
#This function extracts the cds annotation for the gene and exons to be used for getting deseq gene counts 
#and dexseq read counts.
#It returns a list with two objects: 
#1. gene.annotation  <- GRangesList of exons by gene 
#2. exon.annotation <- GRanges of exons
####################################################################################################
ExtractAnnotation <- function(annotation.object.location, genome, num.cores = 12) {

  #Get the ranges corresponding to the coding region
  gr.annotation <- readRDS(annotation.object.location)
  
  if(genome %in% c("hg19", "hg38")){
    vec.chr = c("chrM", "chrX", "chrY", paste("chr", seq(1:22),sep = ""))
  } else if (genome %in% c("mm10")){
  	vec.chr = c("chrMT", "chrX", "chrY", paste("chr", seq(1:19),sep = ""))
  }
  gr.annotation <- gr.annotation[which(seqnames(gr.annotation) %in% vec.chr)]

  grl.gene.annotation <- split(gr.annotation, f = gr.annotation$symbol)
  grl.gene.annotation <- mclapply(grl.gene.annotation, function(gr) {
    values(gr) <- DataFrame(values(gr), exon = paste("exon", seq(1, length(gr)), sep = ""))
    gr
    }, mc.cores = num.cores)
  grl.gene.annotation <- GRangesList(grl.gene.annotation)
  gr.exon.annotation <- unlist(grl.gene.annotation)

  return(list(gene.annotation = grl.gene.annotation, exon.annotation = gr.exon.annotation))

}

######################
ReadRNASeqBam <- function(genome, bam.file, grl.gene.annotation, gr.exon.annotation, 
  paired =  TRUE, num.cores = 1) {
  
  cat("Reading File", bam.file, "\n")
  
  if(genome %in% c("hg19", "hg38")){
    vec.chr = c("chrM", "chrX", "chrY", paste("chr", seq(1:22),sep = ""))
  } else if (genome %in% c("mm10")){
  	vec.chr = c("chrMT", "chrX", "chrY", paste("chr", seq(1:19),sep = ""))
  }

  ls.counts <- mclapply(seq(1,length(vec.chr)), function(x) {
    
    #Specifying range of a particular chromosme
    gr.chr <- GRanges(vec.chr[x], IRanges(start = 1, end = 536870912))
    #Defining the parameters to be read from the bam file
    if(paired){
      parameter <- ScanBamParam(tag=c("NH", "IH"), which = gr.chr, 
        flag = scanBamFlag (isProperPair = TRUE), what=c("flag", "mrnm", "mpos"))      
    } else {
      parameter <- ScanBamParam(tag=c("NH", "IH"), which = gr.chr, 
        what=c("flag", "mrnm", "mpos"))
    }

    ga.alignment <- readGAlignments(bam.file, param = parameter, use.names = TRUE)
    
    #finding unique hits
    if(!all(is.na(values(ga.alignment)$IH))){
      vec.unique.hit <- values(ga.alignment)$IH == 1L
    }

    if(!all(is.na(values(ga.alignment)$NH))){
      vec.unique.hit <- values(ga.alignment)$NH == 1L
    }

    #Get all the reads with unique alignment and then pair it
    ga.unique.alignment <- ga.alignment[which(vec.unique.hit)]

    if(paired){
      gap.alignment <- makeGAlignmentPairs(ga.unique.alignment, use.names = T, use.mcols = T)
      if(length(gap.alignment) == 0){
        stop("No paired reads")
      }
    } else {
      gap.alignment <- ga.unique.alignment
    }

    total.reads <- length(gap.alignment)
    
    #get the counts for overlapping reads
    vec.gene.counts <- unname(assays(summarizeOverlaps(features = grl.gene.annotation, 
      reads = gap.alignment, ignore.strand = T, mode="IntersectionNotEmpty"))$counts)
    
    vec.exon.counts <- unname(assays(summarizeOverlaps(features = gr.exon.annotation, 
      reads = gap.alignment, ignore.strand = T, mode="IntersectionNotEmpty"))$counts)

    rm(list = c("ga.unique.alignment", "gap.alignment"))
    gc()
    return(list(gene = vec.gene.counts, exon = vec.exon.counts, total = total.reads))
  }, mc.cores = num.cores, mc.cleanup = T)
  
  #Sums up the counts across all the chromosomes
  vec.gene.counts <- rowSums(do.call(cbind, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$gene)))
  vec.exon.counts <- rowSums(do.call(cbind, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$exon)))
  vec.total.counts <- sum(do.call(c, lapply(seq(1, length(ls.counts)), function(x) ls.counts[[x]]$total)))

  return(list(gene.counts = vec.gene.counts, exon.counts = vec.exon.counts, 
    library.size = vec.total.counts))
}

####################################################################################################
#This function takes the data table, reads the bam file and gets the corresponsing deseq and dexseq 
#read counts
#It saves two objects: 
#1. gene.annotation  <- GRangesList of exons by gene 
#2. exon.annotation <- GRanges of exons
####################################################################################################

runCDSCounts <- function(project.dir, dt.location, count.obj.dir = NULL) {

  # dt.location <- datainfo.location_uniq

    # --- Initial Setup ---
    if(is.null(count.obj.dir)){
      count.obj.dir <- file.path(project.dir, "rnaseq_counts")
    }
    
    # Create counts directory if needed
    if(!dir.exists(count.obj.dir)){
      dir.create(count.obj.dir, recursive = TRUE)
    }
    
    # --- Data Validation ---
    df.sample.info <- read.delim(dt.location, sep = "\t", header = TRUE)
  
    
    # Ensure unique sample names
    if(any(duplicated(df.sample.info$NAME))) {
      stop("Duplicate entries found in df.sample.info$NAME - all sample names must be unique")
    }
    
    sample.names <- as.character(df.sample.info$NAME)
    
    # --- Sample Processing ---
    results <- sapply(sample.names, function(sample.name) {

      # sample.name <- "MMRF_1038_1_BM"


              # --- Fastq File Handling ---
              fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
             
              
              fastq.files <- unlist(strsplit(fastq.file.location, split = "::"))
              paired <- length(fastq.files) == 2
              
              # Validate fastq files
              # if(!all(file.exists(fastq.files))) {
              #   warning(paste("Missing fastq files for sample:", sample.name))
              #   return(NULL)
              # }
              
              # --- BAM File Handling ---
              bam.file.dir <-as.character(subset(df.sample.info, NAME == sample.name)$FILE_PATH)
            
              
              all.bam.files <- list.files(path = bam.file.dir, pattern = "*.bam$")
              bam.file.logical <- grepl(sample.name, all.bam.files, fixed = TRUE)
              
              if(all(!bam.file.logical)) {
                warning(paste("No BAM files found for sample:", sample.name))

              }
              
              # --- Job Submission ---
              genome.name <- as.character(subset(df.sample.info, NAME == sample.name)$GENOME)
              
              bam.file.location <- file.path(bam.file.dir, all.bam.files[bam.file.logical])
              
              sbatch.cmd <- sprintf(
                paste0("sbatch --export=genome_name=%s,sample_name=%s,bam_file=%s,save_location=%s,paired=%s %s/2_gene_preprocessing/1_rnaseq_pipeline/run_gene_counts.sh"),
                genome.name, sample.name, bam.file.location, count.obj.dir, paired, project.dir
              )

              exit_code <- system(sbatch.cmd)
              if(exit_code != 0) {
                stop(paste("SLURM submission failed for sample:", sample.name))
              }


      })
    
}
