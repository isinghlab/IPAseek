
library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)
library(data.table)



project.dir<-"/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/1_intron_preprocessing/3_filtering_gobj"

setwd(project.dir)

rn_hg38 <- readRDS("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/1_intron_preprocessing/1_flatten_genome/hg38_annotated_numbered_cds.rds")

rn_introns <- readRDS("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/1_intron_preprocessing/2_annotation_object/numbered_hg38_introns_cds.rds")



##############################################################################
                          #Filter Small Introns (<1KB)
##############################################################################

rn_introns_filt_small = rn_introns[which(width(rn_introns) >= 500 & width(rn_introns) <= 150000 )]

saveRDS(rn_introns_filt_small,"rn_introns_filt_small.rds")

##############################################################################
                    # Keep introns in protein coding genes
##############################################################################

rn_cds <- rn_hg38[which(rn_hg38$exon.anno %like% c("cds"))]

rn_cds_id <- unique(rn_cds$entrez.id)

rn_introns_cds <- rn_introns_filt_small[rn_introns_filt_small$entrez.id %in% rn_cds_id]

#rn_introns_cds <- rn_introns[rn_introns$entrez.id %in% rn_cds_id]

write.csv(as.data.frame(rn_introns_cds),"rn_introns_cds.csv")
saveRDS(rn_introns_cds, file = "rn_introns_cds.RDS")


##############################################################################
                             #  Filter sno/miRNA
##############################################################################

sno_miRNA <- read.table(paste0(project.dir,"/filter_files/sno_miRNA"),sep="\t",skip=1)

sno_miRNA_rng<-sno_miRNA[2:4]
names(sno_miRNA_rng)<-c("chr","start","end")
gr_sno_miRNA<-makeGRangesFromDataFrame(sno_miRNA_rng,keep.extra.columns=TRUE)
hits_miRNA = findOverlaps(gr_sno_miRNA,rn_introns_cds,type="within",ignore.strand=TRUE)

index_miRNA<-subjectHits(hits_miRNA)
filt_miRNA_introns<-rn_introns[index_miRNA]
filt_miRNA_introns<-unique(filt_miRNA_introns)
filt_miRNA_introns_id<-names(filt_miRNA_introns)
write.csv(filt_miRNA_introns_id,"filt_miRNA.csv")

##############################################################################
                         # Filter Blacklisted regions
##############################################################################

blacklist<- read.table(paste0(project.dir,"/filter_files/hg38-blacklist"),sep="\t",skip=1)
names(blacklist)<-c("chr","start","end")
gr_blacklist<-makeGRangesFromDataFrame(blacklist ,keep.extra.columns=TRUE)
hits_blacklist = findOverlaps(gr_blacklist,rn_introns_cds,type="within",ignore.strand=TRUE)

index_blacklist<-subjectHits(hits_blacklist)
filt_blacklist_introns<-rn_introns[index_blacklist]
filt_blacklist_introns<-unique(filt_blacklist_introns)
filt_blacklist_introns_id<-names(filt_blacklist_introns)
write.csv(filt_blacklist_introns_id,"filt_blacklisted.csv")

##############################################################################
                          #Filter Retrotransposons
##############################################################################

retrotrans<- read.table(paste0(project.dir,"/filter_files/retroGenes_V9"),sep="\t",skip=1)
drop.cols <- c('V2','V3','V6','V7','V8','V9')
retrotrans <- retrotrans%>% 
dplyr::select(- one_of(drop.cols))
names(retrotrans)<-c("chr","start","end")

retrotrans<-makeGRangesFromDataFrame(retrotrans,keep.extra.columns=TRUE)
hits_retrotrans = findOverlaps(retrotrans,rn_introns_cds,type="within",ignore.strand=TRUE)

index_retrotrans<-subjectHits(hits_retrotrans)
filt_retrotrans_introns<-rn_introns[index_retrotrans]
filt_retrotrans_introns<-unique(filt_retrotrans_introns)
filt_retrotrans_introns_id<-names(filt_retrotrans_introns)
write.csv(filt_retrotrans_introns_id,"filt_retrotransposons.csv")


##############################################################################
                          # Filter the regions
##############################################################################

filt_all_out_id<-unique(c(filt_miRNA_introns_id,filt_blacklist_introns_id,filt_retrotrans_introns_id))

all_introns_id<-names(rn_introns_cds)

selc_intron_ids<-setdiff(all_introns_id,filt_all_out_id)
rn_introns_filt<- rn_introns_cds[ which(names(rn_introns_cds) %in% selc_intron_ids)]

saveRDS(rn_introns_filt,"rnhg38_filtered_introns_cds.rds")

# hg38 <- readRDS("rnhg38_filtered_introns_cds.rds")

# rn_hg38 <- rn_hg38[rn_hg38$entrez.id %in% 3953,]
# rn_introns



