
options(warn = -1)
#########################################################################################
suppressPackageStartupMessages({
  library(yaml)
  library(GenomicFeatures)
  library(Rsamtools)
  library(data.table)
  library(Gviz)
  library(rtracklayer)
  library(GenomicAlignments)
})



##This function takes the location(path) of the bam file, library.size and 
##creates a data track for the bam file
##file.location <- path of the bam file
##color.scheme <- color of the track
##selection <- the granges for which the data track has to be created
generateDataTrackRnaSeq <- function(file.location, color.scheme, selection, 
  vec.library.size, user.ylim){
  
  track.name <- names(file.location)
  gr.data <- myImportFunRnaSeq(file.location, selection)
  gr.data$score <- (gr.data$score)*10^6/vec.library.size[track.name]
  if(is.null(user.ylim)){
    dtrack <- DataTrack(gr.data, col.histogram = color.scheme, fill = color.scheme, 
      name = track.name, genome = "hg19", type = "hist", window = -1, windowSize = 1, 
      baseline = 0, col.baseline= "black", cex.title = 0.4)
  } else {
    dtrack <- DataTrack(gr.data, col.histogram = color.scheme, fill = color.scheme, 
      name = track.name, genome = "hg19", type = "hist", window = -1, windowSize = 1, 
      baseline = 0, col.baseline= "black", cex.title = 0.4, ylim = c(0, user.ylim))
  }
  displayPars(dtrack)$frame <- FALSE
  return(dtrack)
  
}


##This function reads the bam file for the specified GRanges
##It finds the coverage for every nucleotide for the give GRange
##It takes care of the uniqueness of the reads
##file <- bam file
##gr <- GRanges that has to be read
##returns the Granges object with an extra columns: score that gives
##the coverage of the GRanges
myImportFunRnaSeq <- function(file, gr){

  strand(gr) <- "*"
  param <- ScanBamParam(what=c("pos", "qwidth"), which=gr, 
    flag=scanBamFlag(isUnmappedQuery=FALSE), tag = c("NH", "X0"))
  ##Remove the multimapping reads and the antisense reads
  x <- readGAlignments(file, param=param)
  if(!all(is.na(mcols(x)$NH))){
    x <- x[which(mcols(x)$NH == 1)]
    } else if(!all(is.na(mcols(x)$X0))){
      x <- x[which(mcols(x)$X0 == 1)]
    }

  x <- unlist(as(x, "GRangesList"))
  cov <- coverage(IRanges(start(x), width=width(x)))
  if(length(cov)==0){
    gr.data <- gr
    mcols(gr.data) <- DataFrame(score=0)
    }else{
      gr.data <- GRanges(seqnames=seqnames(gr)[1L], 
        ranges=IRanges(start=start(cov), end=end(cov)), strand="*", 
        score=runValue(cov))
      }

  return(gr.data)

}



plotGeneModel <- function(id, 
  vec.library.size = NULL, gr.gene.info, 
  from = NULL, to = NULL,
  plot.dir, 
  plot.rnaseq = F, rna.seq.bams = NULL, highlight.start = NULL, highlight.end = NULL,
  extension = 1000, user.ylim = NULL){

  if(id %in% df.gene.info$entrez.id){

    gr.gene.model <- gr.gene.info[which(gr.gene.info$entrez.id == id)]
    gene.chr <- as.character(seqnames(gr.gene.model)[1L])
    gene.symbol <- gr.gene.model$symbol[1L]
    print(gene.symbol)
    gene.strand <- as.character(strand(gr.gene.model)[1L])
    df.gene.model <- as(gr.gene.model, "data.frame")


    title <- paste(gene.chr, "(", gene.strand , ") " ,gene.symbol, sep = "")

    gtrack <- GenomeAxisTrack()
    grtrack <- GeneRegionTrack(df.gene.model, showId = TRUE, genome = "hg19", collapse = FALSE, 
      name = gene.symbol, chromosome = gene.chr)
    displayPars(grtrack)$fill <- list(fill = "#0000B2")
    displayPars(grtrack)$col <- "#0000B2"
    displayPars(grtrack)$fontcolor.group <- "black"
 
    if(gene.strand == "+"){
      right.extension <- extension
      left.extension <- extension
      } else {
        right.extension <- extension
        left.extension <- extension
      }

    if(is.null(from)){
      from = min(start(gr.gene.model))
    }
    
    if(is.null(to)){
      to = max(end(gr.gene.model))
    }

    from <- from - left.extension
    to <- to + right.extension

    rnaseq.datatracks <- NULL
    gr.gene.rnaseq <- GRanges(seqnames = gene.chr, ranges = IRanges(start = from, end = to),
      strand = gene.strand)
    if(plot.rnaseq && !is.null(rna.seq.bams) & !is.null(vec.library.size)){
      num.files <- length(rna.seq.bams)
      colors <- c("#e41a1c", "#e41a1c", "#377eb8", "#377eb8", "#4daf4a", "#4daf4a",
        "#984ea3", "#984ea3")
      #colors <- rainbow(num.files)
      rnaseq.datatracks <- c()
      for(i in 1:num.files){
        rnaseq.track <- generateDataTrackRnaSeq(file.location = rna.seq.bams[i], color.scheme = colors[i],
          selection = gr.gene.rnaseq, vec.library.size, user.ylim)
        rnaseq.datatracks <- c(rnaseq.datatracks, rnaseq.track)
      }
    }

    ##Creating a highlight track
    real.ht.track <- NULL
    if(!is.null(highlight.start) && !is.null(highlight.end)){
      ht.range <- GRanges(seqnames = gene.chr, ranges = IRanges(start = highlight.start, end = highlight.end))
      real.ht.track <- HighlightTrack(trackList = data.tracks, range = ht.range, genome = "hg19")
      displayPars(real.ht.track)$fill <- "#efefef"
      displayPars(real.ht.track)$col <- "#efefef"
    }
   
    plotTracks(c(rnaseq.datatracks, grtrack, gtrack), from = from, to = to, 
      col.title = "black", collapse = FALSE, col.axis = "black", 
      background.title = "transparent", main = title, cex.main = 1)


    } else {
    sprintf("Gene model for entrez id %s is missing", id)
  }

}
