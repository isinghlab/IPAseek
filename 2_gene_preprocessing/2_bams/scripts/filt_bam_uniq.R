
filt_bam<-function(dt.location,wd){


slurm.files.dir <- file.path(wd, "slurm_submission")
logs.files.dir <- file.path(wd, "logs")


if(!dir.exists(slurm.files.dir )){
        dir.create(slurm.files.dir )
   }




data.input <- read.delim(dt.location, sep="\t", header=T)
sampleNames <- as.character(data.input$NAME)


#Testing samples
#sampleNames<-c("RNA-seq_CLL7_2")


sapply(sampleNames, function(sample){

#store the path for bam file
bam.loc1 <- as.character(subset(data.input, NAME==sample)$FILE_PATH)
bam.loc2 <- as.character(subset(data.input, NAME==sample)$UNIQUE_ID)
    
bam.loc <- paste0(bam.loc1,"/",bam.loc2,".bam")



#store the bam file directory
bam.dir <- dirname(bam.loc)
#in case bam file is missing show an error message
if(!file.exists(bam.loc)){
    print("ERROR! Bam file does not exist!")
}





paths <- data.frame(bam.loc)
paths$bam.dir <- bam.dir
paths$sample <- sample

filtBamUniq.loc<-paste0(slurm.files.dir,"/filtBamUniq.Rdata")
save("filtBamUniq", file=filtBamUniq.loc)
    
if(!dir.exists(paste0(slurm.files.dir,"/", sample))){
          dir.create(paste0(slurm.files.dir,"/", sample))
     }
if(!dir.exists(logs.files.dir)){
          dir.create(logs.files.dir)
    }




sample.dir <- file.path(paste0(slurm.files.dir,"/", sample, "/"))
script.name <- file.path(paste0(sample.dir,'/', "_filtBamUniq.R"))
save(paths, file=paste0(sample.dir,sample,"_","paths.Rdata"))
sink(file=script.name)



cat(paste0("\nload (\'",filtBamUniq.loc,"')"))
#load the previously saved paths object
cat(paste0("\nload (\'",sample.dir,sample,"_", "paths.Rdata\')"))
#call the function
cat(paste0("\nfiltBamUniq(paths)"))
#close
sink()

#create the bash script
      bash.file.location <- file.path(paste0(sample.dir,'submit.sh'))
      slurm.jobname <- sprintf("#SBATCH --job-name=%s_filt_bam", sample)
      slurm.time <- sprintf("#SBATCH --time=01:30:00")
      slurm.mem <- sprintf("#SBATCH --mem=12G")
      slurm.tasks <- sprintf("#SBATCH --ntasks=8")
      slurm.output <- paste0("#SBATCH --output=",logs.files.dir,"/",sample,"_filt_bam")
      sbatch.line<- sprintf("Rscript --vanilla %s", script.name)
      file.conn <- file(bash.file.location)
      writeLines(c("#!/bin/bash", slurm.jobname, slurm.time, slurm.mem, slurm.tasks, slurm.output, sbatch.line), file.conn)
      close(file.conn)
      system(paste0("sbatch ", bash.file.location))

   
 print(paste0("Running ... ", "job for ", sample, " submitted."))



})

}


filtBamUniq<-function(paths){

start_time <- Sys.time()


print("started filtering bam files")

bam.loc<-paths$bam.loc
bam.dir<-paths$bam.dir
sample<-paths$sample


uniq.bam.dir <- file.path(bam.dir, "uniq_bams")


if(!dir.exists(uniq.bam.dir)){
        dir.create(uniq.bam.dir)
   }

#uniquely mapped reads bam file loc
uniq.bam.loc <- paste0(uniq.bam.dir,"/",sample,"_uniq.bam")
uniq.bam.loc.index<-paste0(uniq.bam.dir,"/",sample,"_uniq.bam.bai")



#if the uniquely mapped bam file does not exist
if(!file.exists(uniq.bam.loc)){
print("filter")
#filter for uniquely mapped reads from bam file
system(paste0("samtools view -h ", bam.loc," | grep -P \"(NH:i:1|^@)\" | samtools view -q 255 -b > ", uniq.bam.dir, "/", sample, "_uniq.bam"))

print(paste0("unique bam for ",sample," is created...... wait generating index file"))
}

if(!file.exists(uniq.bam.loc.index)){
#create the index file for uniquely mapped reads
system(paste0("samtools index ", uniq.bam.dir, "/",sample, "_uniq.bam ",uniq.bam.dir, "/", sample, "_uniq.bam.bai"))

print(paste0("index file generated for ",sample))
}



end_time <- Sys.time()
exec.time<-end_time - start_time
print(paste0("execution time:",exec.time))

}



