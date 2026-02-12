
library(parallel)

calc_atlas <- function(data.input){

   #data.input <- data.input.split
   data.input_name <-data.input$atlas[1]
   wd<- data.input$wd[1]
   atlas_name<- data.input$atlas_name[1]

   source(paste0(wd,"/3_ipa_run/scripts/7_calc_cds_overlap.R"))

    
   groups_df <- aggregate(NAME ~ CELL_TYPE, data.input[c("CELL_TYPE", "NAME")], function(x) paste(x, collapse = ","))


   #Get the sample names
   sampleNames <- as.character(data.input$NAME)
   source_name <- as.character(data.input$ATLAS_NAME)
# 
#    print(paste("processing for :", sampleNames))
# 	 
# 	# cpt_all<-do.call(rbind, lapply(list.files(path = paste0(wd,"/pelt/results/",source_name,"/",sampleNames,"/exon_exprs_results"),full.names = TRUE, pattern = "*.csv") , read.csv))
# 	#  
# 	cpt_all <- do.call(
# 	  rbind,
# 	  lapply(
# 	    # List all te_expression*.csv files
# 	    list.files(
# 	      path = paste0(wd, "/pelt/results/", source_name, "/", sampleNames, "/exon_exprs_results/"),
# 	      full.names = TRUE,
# 	      pattern = "^te_expression.*\\.csv$"
# 	    ) |> 
# 	      # Exclude files ending with _unfiltered.csv
# 	      (\(x) x[!grepl("_unfiltered\\.csv$", x)] )(),
# 	    read.csv
# 	  )
# 	)
# 	
# 	
# 	  cpt_all <- unique(cpt_all)
#     # cpt_all <- cpt_all[-1]
#     cpt_all <- cpt_all[c("seqnames", "start", "end", "strand", "entrez.id", "id", "ipa_sel", "source")]
# 
#     cpt_all$source_group <- data.input$CONDITION[match(cpt_all$source,data.input$NAME)]
# 
#     # make granges object for the terminal exon
#     te_gr <- makeGRangesFromDataFrame(cpt_all, keep.extra.columns=TRUE)
# 
#     # te_gr <- head(te_gr)
# 
#     te_gr$entrez.id <- as.character(te_gr$entrez.id)
#     te_gr$split_col <- paste(te_gr$id, te_gr$ipa_sel, sep=" ")
# 
#     # split the granges object for each IPA 
#     te_gr_split<-split(te_gr, f = c(te_gr$split_col))
# 
#     #te_gr_split <- te_gr_split[1:50]
# 
#     # reduce the IPA within 100bp into one
#     print(paste("Making Atlas"))
# 
#     te_median <- lapply(te_gr_split, function(gr){
#       
#       gr_orig <- gr                        # Store the original group
#       gr <- reduce(gr, min.gapwidth=100L)  # Merge regions within 100bp
#       
#       gr <- annotateGRfromGR(gr, gr_orig)  # Annotate merged region with original metadata
#       
#       # For each merged region, set the terminal exon boundary to the median of the original ends
#       gr_split <- split(gr)
#       gr_median <- lapply(gr_split, function(gr_splitted){
#         
#         # For positive strand, use the median of the ends
#         if(as.character(strand(gr_splitted)) == "+" ){
#           gr_within <- subsetByOverlaps(gr_orig, gr_splitted)
#           ends <- c(end(gr_within))
#           ends_median <- median(ends)
#           end(gr_splitted) <- ends_median
#         } else { # For negative strand, use the median of the starts
#           gr_within <- subsetByOverlaps(gr_orig, gr_splitted)
#           ends <- c(start(gr_within))
#           ends_median <- median(ends)
#           start(gr_splitted) <- ends_median
#         }
#         return(gr_splitted)
#       })
#       gr_median <- unlist(gr_median)
#       gr_median <- GRangesList(gr_median)
#       gr_median <- unlist(gr_median)
#       print(gr_median) # Print progress for each group
#     })
#       
#     
#     print("calculating retained cds")
#     te_median <- lapply(te_median,calc_cds_overlap)
# 
#     te_median<-GRangesList(te_median)
#     te_median<-unlist(te_median)
# 
# 
#     te_median$no_of_sources <- sapply(strsplit(te_median$source,','), uniqueN)
# 
#     te_median$percentage_of_sources <- (te_median$no_of_sources/length(sampleNames))*100
        res.dir <- file.path(wd,"pelt", "results", atlas_name)
# 
#     if(!dir.exists(res.dir )){
#         dir.create(res.dir )
#       }


    te_median <- readRDS(paste0("/scratch/user/richa.rashmi.1202/ipa/IPAseek_pipeline/pelt/results/combined/combined_full_ipa_atlas.RDS"))

    te_median_confident <- calc_confidence(te_median,groups_df)


    print(paste("Saving Atlas"))
    print(res.dir)
    print(atlas_name)
    print(data.input_name)

    # saveRDS(te_median, paste0(res.dir, "/",atlas_name,"_",data.input_name, "_ipa_atlas.RDS"))
    # write.csv(as.data.frame(te_median), paste0(res.dir, "/",atlas_name,"_",data.input_name,"_ipa_atlas.csv"))

    saveRDS(te_median_confident, paste0(res.dir, "/",atlas_name,"_",data.input_name, "_ipa_atlas_conf.RDS"))
    write.csv(as.data.frame(te_median_confident), paste0(res.dir, "/",atlas_name,"_",data.input_name,"_ipa_atlas_conf.csv"))

   

}