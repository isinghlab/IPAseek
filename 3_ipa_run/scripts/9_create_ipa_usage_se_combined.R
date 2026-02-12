

library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)
library(tidyr)
library(data.table)
library(gtools)




ipa_usage_se<-function(input.data.path,wd, atlas_name){


# input.data.path <- "/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/data_tables/combined.txt"
# wd <- "/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline"
# atlas_name <- "combined"

## Get the input table with sample info 
data.input <- read.delim(input.data.path, sep="\t", header=T)

source_name <- unique(data.input$SOURCE_NAME)

## Get the IPA atlas

ipa_atlas <-  readRDS(paste0(wd,"/pelt/results/",atlas_name,"/", atlas_name,"_full_ipa_atlas_conf.RDS"))

ipa_atlas  <- split(ipa_atlas, as.factor(ipa_atlas))
ipa_atlas<-GRangesList(ipa_atlas)
ipa_atlas<-unlist(ipa_atlas)


## Get the samples included in the atlas

  sampleNames <- as.character(data.input$NAME)

## Define assays for the IPA usage

  # utr3_ipa_usage <- data.frame()
  cds_ipa_usage <- data.frame()
  cdsex100_ipa_usage <- data.frame()
  cdsex200_ipa_usage <- data.frame()
  cdsex300_ipa_usage <- data.frame()
  cdsex500_ipa_usage <- data.frame()
  # cpe_ipa_usage <- data.frame()

  ipa_reads <- data.frame()
  cds_reads <- data.frame()
  cdsex100_reads <- data.frame()
  cdsex200_reads <- data.frame()
  cdsex300_reads <- data.frame()
  cdsex500_reads <- data.frame()
  # utr3_reads <- data.frame()
  # cpe_reads <- data.frame()

  ipa_tpm <- data.frame()
  cds_tpm <- data.frame()
  cdsex100_tpm <- data.frame()
  cdsex200_tpm <- data.frame()
  cdsex300_tpm <- data.frame()  
  cdsex500_tpm <- data.frame() 
  # utr3_tpm <- data.frame()
  # cpe_tpm <- data.frame()

## get information from different samples and merge them into the assays

   for(sample in sampleNames) {

      # sample <-"ENCFF795CTQ"
    source_name <- data.input[data.input$NAME == sample,"ATLAS_NAME"]

    x <- read.csv(paste0(wd,"/pelt/results/",atlas_name,"/",source_name,"_", sample,"_ipa_usage_atlas.csv"))
    x[is.na(x)] = 0
    print(head(x))

    

    #x <- x[x$percentage_of_sources >25,]

    # x_utr3_ipa_usage <- x[c("X", "utr3_ipa_usage")]
    # setnames(x_utr3_ipa_usage, old = c("utr3_ipa_usage") , new =sample )

    x_cds_ipa_usage <- x[c("X", "cds_ipa_usage")]
    setnames(x_cds_ipa_usage, old = c("cds_ipa_usage") , new =sample )

    x_cdsex100_ipa_usage <- x[c("X", "cdsex100_ipa_usage")]
    setnames(x_cdsex100_ipa_usage, old = c("cdsex100_ipa_usage") , new =sample )

    x_cdsex200_ipa_usage <- x[c("X", "cdsex200_ipa_usage")]
    setnames(x_cdsex200_ipa_usage, old = c("cdsex200_ipa_usage") , new =sample )

    x_cdsex300_ipa_usage <- x[c("X", "cdsex300_ipa_usage")]
    setnames(x_cdsex300_ipa_usage, old = c("cdsex300_ipa_usage") , new =sample )

    x_cdsex500_ipa_usage <- x[c("X", "cdsex500_ipa_usage")]
    setnames(x_cdsex500_ipa_usage, old = c("cdsex500_ipa_usage") , new =sample )

    # x_cpe_ipa_usage <- x[c("X", "cpe_ipa_usage")]
    # setnames(x_cpe_ipa_usage, old = c("cpe_ipa_usage") , new =sample )


    x_ipa_reads <- x[c("X", "ipa_reads")]
    setnames(x_ipa_reads, old = c("ipa_reads") , new =sample )

    x_cds_reads <- x[c("X", "cds_reads")]
    setnames(x_cds_reads, old = c("cds_reads") , new =sample )

    x_cdsex100_reads <- x[c("X", "cdsex100_reads")]
    setnames(x_cdsex100_reads, old = c("cdsex100_reads") , new =sample )

    x_cdsex200_reads <- x[c("X", "cdsex200_reads")]
    setnames(x_cdsex200_reads, old = c("cdsex200_reads") , new =sample )

    x_cdsex300_reads <- x[c("X", "cdsex300_reads")]
    setnames(x_cdsex300_reads, old = c("cdsex300_reads") , new =sample )

    x_cdsex500_reads <- x[c("X", "cdsex500_reads")]
    setnames(x_cdsex500_reads, old = c("cdsex500_reads") , new =sample )

    # x_utr3_reads <- x[c("X", "utr3_reads")]
    # setnames(x_utr3_reads, old = c("utr3_reads") , new =sample )

    # x_cpe_reads <- x[c("X", "cpe_reads")]
    # setnames(x_cpe_reads, old = c("cpe_reads") , new =sample )
    

    x_ipa_tpm <- x[c("X", "ipa_tpm")]
    setnames(x_ipa_tpm, old = c("ipa_tpm") , new =sample )

    x_cds_tpm <- x[c("X", "cds_tpm")]
  	setnames(x_cds_tpm, old = c("cds_tpm") , new =sample )

    x_cdsex100_tpm <- x[c("X", "cdsex100_tpm")]
    setnames(x_cdsex100_tpm, old = c("cdsex100_tpm") , new =sample )

    x_cdsex200_tpm <- x[c("X", "cdsex200_tpm")]
    setnames(x_cdsex200_tpm, old = c("cdsex200_tpm") , new =sample )

    x_cdsex300_tpm <- x[c("X", "cdsex300_tpm")]
    setnames(x_cdsex300_tpm, old = c("cdsex300_tpm") , new =sample )

    x_cdsex500_tpm <- x[c("X", "cdsex500_tpm")]
    setnames(x_cdsex500_tpm, old = c("cdsex500_tpm") , new =sample )

    # x_utr3_tpm <- x[c("X", "utr3_tpm")]
    # setnames(x_utr3_tpm, old = c("utr3_tpm") , new =sample )

    # x_cpe_tpm <- x[c("X", "cpe_tpm")]
    # setnames(x_cpe_tpm, old = c("cpe_tpm") , new =sample )




   if (nrow(cds_ipa_usage) == 0) {

    # utr3_ipa_usage <- x_utr3_ipa_usage
    cds_ipa_usage <- x_cds_ipa_usage
    cdsex100_ipa_usage <- x_cdsex100_ipa_usage
    cdsex200_ipa_usage <- x_cdsex200_ipa_usage
    cdsex300_ipa_usage <- x_cdsex300_ipa_usage
    cdsex500_ipa_usage <- x_cdsex500_ipa_usage
    # cpe_ipa_usage <- x_cpe_ipa_usage

    ipa_reads <- x_ipa_reads
	  cds_reads <- x_cds_reads
    cdsex100_reads <- x_cdsex100_reads
    cdsex200_reads <- x_cdsex200_reads
    cdsex300_reads <- x_cdsex300_reads
    cdsex500_reads <- x_cdsex500_reads
    # utr3_reads <- x_utr3_reads
    # cpe_reads <- x_cpe_reads

	  ipa_tpm <- x_ipa_tpm
	  cds_tpm <- x_cds_tpm
    cdsex100_tpm <- x_cdsex100_tpm
    cdsex200_tpm <- x_cdsex200_tpm
    cdsex300_tpm <- x_cdsex300_tpm
    cdsex500_tpm <- x_cdsex500_tpm
    # utr3_tpm <- x_utr3_tpm
    # cpe_tpm <- x_cpe_tpm
    
    } else {

      # utr3_ipa_usage <- merge(utr3_ipa_usage, x_utr3_ipa_usage)
      cds_ipa_usage <- merge(cds_ipa_usage, x_cds_ipa_usage)
      cdsex100_ipa_usage <- merge(cdsex100_ipa_usage, x_cdsex100_ipa_usage)
      cdsex200_ipa_usage <- merge(cdsex200_ipa_usage, x_cdsex200_ipa_usage)
      cdsex300_ipa_usage <- merge(cdsex300_ipa_usage, x_cdsex300_ipa_usage)
      cdsex500_ipa_usage <- merge(cdsex500_ipa_usage, x_cdsex500_ipa_usage)
      # cpe_ipa_usage <- merge(cpe_ipa_usage, x_cpe_ipa_usage)

      ipa_reads <- merge(ipa_reads, x_ipa_reads)
      cds_reads <- merge(cds_reads, x_cds_reads)
      cdsex100_reads <- merge(cdsex100_reads, x_cdsex100_reads)
      cdsex200_reads <- merge(cdsex200_reads, x_cdsex200_reads)
      cdsex300_reads <- merge(cdsex300_reads, x_cdsex300_reads)
      cdsex500_reads <- merge(cdsex500_reads, x_cdsex500_reads)
      # utr3_reads <- merge(utr3_reads, x_utr3_reads)
      # cpe_reads <- merge(cpe_reads, x_cpe_reads)

	    ipa_tpm <- merge(ipa_tpm, x_ipa_tpm)
      cds_tpm <- merge(cds_tpm, x_cds_tpm)
      cdsex100_tpm <- merge(cdsex100_tpm, x_cdsex100_tpm)
      cdsex200_tpm <- merge(cdsex200_tpm, x_cdsex200_tpm)
      cdsex300_tpm <- merge(cdsex300_tpm, x_cdsex300_tpm)
      cdsex500_tpm <- merge(cdsex500_tpm, x_cdsex500_tpm)
      # utr3_tpm <- merge(utr3_tpm, x_utr3_tpm)
      # cpe_tpm <- merge(cpe_tpm, x_cpe_tpm)

     #print(head(ipa_file))
     #print(head(ipa_row_data))
    }

   }

   # rownames(utr3_ipa_usage) <- utr3_ipa_usage$X
   # utr3_ipa_usage<- utr3_ipa_usage[order(row.names(utr3_ipa_usage)), ]
   # utr3_ipa_usage <- utr3_ipa_usage[-1]

   rownames(cds_ipa_usage) <- cds_ipa_usage$X
   cds_ipa_usage<- cds_ipa_usage[order(row.names(cds_ipa_usage)), ]
   cds_ipa_usage <- cds_ipa_usage[-1]

   rownames(cdsex100_ipa_usage) <- cdsex100_ipa_usage$X
   cdsex100_ipa_usage<- cdsex100_ipa_usage[order(row.names(cdsex100_ipa_usage)), ]
   cdsex100_ipa_usage <- cdsex100_ipa_usage[-1]

   rownames(cdsex200_ipa_usage) <- cdsex200_ipa_usage$X
   cdsex200_ipa_usage<- cdsex200_ipa_usage[order(row.names(cdsex200_ipa_usage)), ]
   cdsex200_ipa_usage <- cdsex200_ipa_usage[-1]

   rownames(cdsex300_ipa_usage) <- cdsex300_ipa_usage$X
   cdsex300_ipa_usage<- cdsex300_ipa_usage[order(row.names(cdsex300_ipa_usage)), ]
   cdsex300_ipa_usage <- cdsex300_ipa_usage[-1]

   rownames(cdsex500_ipa_usage) <- cdsex500_ipa_usage$X
   cdsex500_ipa_usage<- cdsex500_ipa_usage[order(row.names(cdsex500_ipa_usage)), ]
   cdsex500_ipa_usage <- cdsex500_ipa_usage[-1]

   # rownames(cpe_ipa_usage) <- cpe_ipa_usage$X
   # cpe_ipa_usage<- cpe_ipa_usage[order(row.names(cpe_ipa_usage)), ]
   # cpe_ipa_usage <- cpe_ipa_usage[-1]

   rownames(ipa_reads) <- ipa_reads$X
   ipa_reads<- ipa_reads[ order(row.names(ipa_reads)), ]
   ipa_reads <- ipa_reads[-1]

   rownames(cds_reads) <- cds_reads$X
   cds_reads<- cds_reads[ order(row.names(cds_reads)), ]
   cds_reads <- cds_reads[-1]

   rownames(cdsex100_reads) <- cdsex100_reads$X
   cdsex100_reads<- cdsex100_reads[ order(row.names(cdsex100_reads)), ]
   cdsex100_reads <- cdsex100_reads[-1]

   rownames(cdsex200_reads) <- cdsex200_reads$X
   cdsex200_reads<- cdsex200_reads[ order(row.names(cdsex200_reads)), ]
   cdsex200_reads <- cdsex200_reads[-1]

   rownames(cdsex300_reads) <- cdsex300_reads$X
   cdsex300_reads<- cdsex300_reads[ order(row.names(cdsex300_reads)), ]
   cdsex300_reads <- cdsex300_reads[-1]

   rownames(cdsex500_reads) <- cdsex500_reads$X
   cdsex500_reads<- cdsex500_reads[ order(row.names(cdsex500_reads)), ]
   cdsex500_reads <- cdsex500_reads[-1]

   # rownames(utr3_reads) <- utr3_reads$X
   # utr3_reads<- utr3_reads[ order(row.names(utr3_reads)), ]
   # utr3_reads <- utr3_reads[-1]

   # rownames(cpe_reads) <- cpe_reads$X
   # cpe_reads<- cpe_reads[ order(row.names(cpe_reads)), ]
   # cpe_reads <- cpe_reads[-1]

   rownames(ipa_tpm) <- ipa_tpm$X
   ipa_tpm<- ipa_tpm[ order(row.names(ipa_tpm)), ]
   ipa_tpm <- ipa_tpm[-1]

   rownames(cds_tpm) <- cds_tpm$X
   cds_tpm<- cds_tpm[ order(row.names(cds_tpm)), ]
   cds_tpm <- cds_tpm[-1]

   rownames(cdsex100_tpm) <- cdsex100_tpm$X
   cdsex100_tpm<- cdsex100_tpm[ order(row.names(cdsex100_tpm)), ]
   cdsex100_tpm <- cdsex100_tpm[-1]

   rownames(cdsex200_tpm) <- cdsex200_tpm$X
   cdsex200_tpm<- cdsex200_tpm[ order(row.names(cdsex200_tpm)), ]
   cdsex200_tpm <- cdsex200_tpm[-1]

   rownames(cdsex300_tpm) <- cdsex300_tpm$X
   cdsex300_tpm<- cdsex300_tpm[ order(row.names(cdsex300_tpm)), ]
   cdsex300_tpm <- cdsex300_tpm[-1]

   rownames(cdsex500_tpm) <- cdsex500_tpm$X
   cdsex500_tpm<- cdsex500_tpm[ order(row.names(cdsex500_tpm)), ]
   cdsex500_tpm <- cdsex500_tpm[-1]

   # rownames(utr3_tpm) <- utr3_tpm$X
   # utr3_tpm<- utr3_tpm[ order(row.names(utr3_tpm)), ]
   # utr3_tpm <- utr3_tpm[-1]

   # rownames(cpe_tpm) <- cpe_tpm$X
   # cpe_tpm<- cpe_tpm[ order(row.names(cpe_tpm)), ]
   # cpe_tpm <- cpe_tpm[-1]


   ## create the rowRanges object

   ipa_atlas<- ipa_atlas[ order(names(ipa_atlas)), ]
   rowRanges_ipa <- ipa_atlas[names(ipa_atlas) %in% rownames(cds_ipa_usage),]

   ## create the colData object
  	
   colData_ipa <- data.input
   rownames(colData_ipa) <- data.input$NAME

   ## create the SummarizedExperiment object
   ipa.se <- SummarizedExperiment(assays=list( cds_ipa_usage=cds_ipa_usage,
                                              cdsex100_ipa_usage=cdsex100_ipa_usage,
                                              cdsex200_ipa_usage=cdsex200_ipa_usage,
                                              cdsex300_ipa_usage=cdsex300_ipa_usage,
                                              cdsex500_ipa_usage=cdsex500_ipa_usage,
                                              # utr3_ipa_usage=utr3_ipa_usage,
                                              # cpe_ipa_usage=cpe_ipa_usage,
                                              ipa_reads=ipa_reads,
                                              cds_reads=cds_reads,
                                              cdsex100_reads=cdsex100_reads,
                                              cdsex200_reads=cdsex200_reads,
                                              cdsex300_reads=cdsex300_reads,
                                              cdsex500_reads=cdsex500_reads,
                                              # utr3_reads=utr3_reads,
                                              # cpe_reads=cpe_reads,
                                              ipa_tpm=ipa_tpm,
                                              cds_tpm=cds_tpm,
                                              cdsex_tpm=cdsex100_tpm,
                                              cdsex_tpm=cdsex200_tpm,
                                              cdsex_tpm=cdsex300_tpm,
                                              cdsex_tpm=cdsex500_tpm
                                              # utr3_tpm=utr3_tpm,
                                              # cpe_tpm=cpe_tpm   
                                              ), colData=colData_ipa, rowRanges=rowRanges_ipa)

   ## save the object
   print(paste("Saving SE for : ", sample))
   saveRDS(ipa.se,paste0("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/pelt/results/", atlas_name,"/",atlas_name, "_ipa_usage_se.rds"))
   
   # ipa.se <- readRDS(paste0("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/pelt/results/", atlas_name,"/",atlas_name, "_ipa_usage_se.rds"))
   
   # assay(ipa.se[,ipa.se$CONDITION ==  "no_pma"])

   print("finished")


}


