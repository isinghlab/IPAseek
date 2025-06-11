library(GenomicFeatures)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library(ggplot2)
library(NMF)
library(RColorBrewer)
library(data.table)
library(GenomicRanges)
library(preprocessCore)
library(org.Hs.eg.db)
getExpressedGenes <- function(se.gene, reqd.cond = NULL, rpkm.cutoff = 1, exprsd.in = 0.5){

	coldata.clmns <- colnames(colData(se.gene))
	if(!("condition" %in% coldata.clmns)) stop("condition column absent")

	all.cond <- unique(colData(se.gene)$condition)
	if(!is.null(reqd.cond)){
		reqd.cond <- intersect(all.cond, reqd.cond)
	} else {
		reqd.cond <- all.cond
	}
	
	if(length(reqd.cond) > 0){
		m.logicals <- sapply(reqd.cond, function(idx) {
		  print(idx)
		  col.pos <- which(colData(se.gene)$condition == idx)
		  m.rpkm <- as.matrix(assays(se.gene)$rpkm[ ,col.pos])
		  apply(m.rpkm, 1, function(x) sum(x > rpkm.cutoff) > (exprsd.in * length(x)))
		})
		exprs <- apply(m.logicals, 1, any)
		se.gene.exprs <- se.gene[exprs]	
	}


}

createSE <- function(project.dir, count.files.dir = NULL, data.table.name, design.colmns = NULL){

	datainfo.location <- file.path(project.dir, "data_tables", data.table.name)
	df.design <- read.delim(datainfo.location, sep = "\t", header  = T, stringsAsFactors= F)
	rownames(df.design) <- df.design$NAME

	if(is.null(count.files.dir)){
		count.files.dir <- file.path(project.dir, "rnaseq_counts")
	}
	
	counts.files <- list.files(path = count.files.dir, pattern = "*_gene_counts.rds", 
	    full.names = T)
	expt.names <- gsub(pattern = "_gene_counts.rds", replacement = "", basename(counts.files))

	reqd.counts.files <- counts.files[which(expt.names %in% rownames(df.design))]
	if(length(reqd.counts.files) != length(rownames(df.design))){
		stop("Some experiment count is missing")
	}

	se.gene <- combineSE(reqd.counts.files)
	m.rpkm <- RPKM(se.gene, from.se = F)
	assays(se.gene)$rpkm <- m.rpkm

	temp <- colData(se.gene)
	if(!is.null(design.colmns)){
		df <- df.design[colnames(se.gene), design.colmns, drop = F]
		colData(se.gene) <- DataFrame(temp, df)
	}

	if("CONDITION" %in% design.colmns){
		colData(se.gene)$condition <- colData(se.gene)$CONDITION
		colData(se.gene)$CONDITION <- NULL
	}

	se.gene

}

##############################################################
##This combines the individual summarized experiment objects to 
##create one summarized experiment object
##############################################################
combineSE <- function(rds.files){

  grl.gene.annotation <- rowRanges(readRDS(rds.files[1]))

  ##get all the counts for the samples
  m.count <- {}
  expt.names <- {}
  for(rds.file in rds.files){
    se.expt <- readRDS(rds.file)
    m.count <- cbind(m.count, assays(se.expt)$count)
    expt.names <- c(expt.names, colnames(se.expt))
  }
  colnames(m.count) <- expt.names

  se.gene <- SummarizedExperiment(rowData = grl.gene.annotation, 
    assays = SimpleList(count = m.count))

  ##If library size is given then add that
  coldata.cols <- names(colData(readRDS(rds.files[1])))
  if("library.size" %in% coldata.cols){
    lib.size <- {}
    for(rds.file in rds.files){
      se.expt <- readRDS(rds.file)
      lib.size <- c(lib.size, colData(se.expt)$library.size)
    }
    df.lib.size <- DataFrame(library.size = lib.size)
    rownames(df.lib.size) <- expt.names
    colData(se.gene) <- df.lib.size
  }

  se.gene
  
}

getGeneSE <- function(se.gene){
	gr.exon.ranges <- unlist(rowRanges(se.gene))
	names(gr.exon.ranges) <- NULL
	dt.exon.ranges <- as(as(gr.exon.ranges, "data.frame"), "data.table")
	dt.gene.ranges <- dt.exon.ranges[,{
	  start <- min(start)
	  end <- max(end)
	  symbol <- symbol[1L]
	  seqnames <- seqnames[1L]
	  xx <- list(start = start, end = end, symbol = symbol, seqnames= seqnames)
	  }, by = "gene_id"]
	gr.gene <- as(as(dt.gene.ranges, "data.frame"), "GRanges")
	rowRanges(se.gene) <- gr.gene

	se.gene
}

doDiffExpByDESeq <- function(se.gene, all.enames, cmpr.expts = NULL, count.idx = 1L, 
	condition.clmn = "condition") {

	if(!is.null(all.enames) && (all.enames %in% colnames(se.gene))){
		se.gene. <- se.gene[ ,all.enames]
	}

	ref.expt <- cmpr.expts[1L]
	other.expts <- setdiff(cmpr.expts, ref.expt)

	dds <- DESeqDataSetFromMatrix(countData = assays(se.gene)[["count"]],
                              colData = colData(se.gene),
                              design = ~condition)

	df.results <- do.call(cbind, lapply(other.expts, function(expt){
		dds.sub <- dds[ , which(colData(dds)$condition %in% c(ref.expt, expt))]
		dds.sub$condition <- droplevels(dds.sub$condition)
		dds.sub$condition <- relevel(dds.sub$condition, ref = ref.expt)
		dds.sub <- DESeq(dds.sub)

    comparison.done <- resultsNames(dds.sub)[2]
    comparison.done <- gsub(pattern = sprintf("%s_", condition.clmn), replacement = "", x = comparison.done)
    print(comparison.done)

		df.res <- as.data.frame(results(dds.sub))
		comparison.name <- sprintf("%s_vs_%s", expt, ref.expt)
		old.clmn.names <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
		new.clmn.names <- sprintf("%s_%s", comparison.done, old.clmn.names)
		setnames(df.res, old = old.clmn.names, new = new.clmn.names)
		return(df.res)
		}))

	mcols(se.gene) <- DataFrame(cbind(mcols(se.gene), df.results))
	return(se.gene)
}

writeFile <- function(df, file.name){
	write.table(df, file.name, sep = "\t", col.names = F, row.names = F, quote = F)
}

gseaRankedList <- function(df, cmp.name, gene.column = "symbol", file.location = NULL){

	fc.name <- paste(cmp.name, "log2FoldChange", sep = "_")
	pvalue.name <- paste(cmp.name, "pvalue", sep = "_")

	df.diff <- df[ ,c("gene_id", "symbol", fc.name, pvalue.name)]
	df.diff$rank <- -log10(df.diff[[pvalue.name]]) * sign(df.diff[[fc.name]])
	df.diff <- df.diff[which(!is.na(df.diff$rank)), ]

	max.vals <- unique(sort(df.diff$rank, decreasing = T))
	if(is.infinite(max.vals[1L])){
	  max.val <- max.vals[2L]
	  } else {
	    max.val <- max.vals[1L]
	  }
	min.vals <- unique(sort(df.diff$rank))
	if(is.infinite(min.vals[1L])){
	  min.val <- min.vals[2L]
	  } else {
	    min.val <- min.vals[1L]
	  }
	df.diff$rank <- ifelse(df.diff$rank > max.val, max.val, df.diff$rank)
	df.diff$rank <- ifelse(df.diff$rank < min.val, min.val, df.diff$rank)
	
	df.rnkd <- df.diff[ ,c(gene.column, "rank")]
	df.rnkd <- df.rnkd[order(df.rnkd$rank, decreasing = T), ]
	colnames(df.rnkd) <- c("id", "rank")

  if(!is.null(file.location)){
    write.table(df.rnkd, file = file.location, 
      col.names = T, row.names = F, quote = F, sep = "\t")    
  }

  return(df.rnkd)

}

performGO <- function(ids, universe, annotation=NULL){
    gobp <- do.gostats(ids, universe, ontology=c('BP'),
        conditional=FALSE, p.value=0.05, testDirection='over',
        annotation)
    df.gobp <- summary(gobp)
    return(df.gobp)
}

writeReport <- function(diff.genes, bkg.genes, ann = NULL, file.location){
	x <- performGO(diff.genes, bkg.genes, annotation=ann)
	write.table(subset(x, Size <= 400 & Size >= 10), 
    	file = file.location, 
    	col.names = T, row.names = F, quote = F, sep = "\t")
}

generateVolcanoPlot <- function(go.id = NULL, df.diff, comp.name = NULL, 
  pval.plot.cutoff = NULL, FC.plot.cutoff = NULL,
  fc.cutoff  = 2, log10adj.cutoff = 2, label = TRUE, label.fc.cutoff  = 2, 
  point.size = 1, label.size = 4, label.log10adj.cutoff  = 2){

	if(!is.null(go.id)){
    	xx <- as.list(org.Hs.egGO2ALLEGS)
    	go.gene.ids <- xx[[go.id]]
    	df.hs.symbol <- toTable(org.Hs.egSYMBOL)
	}

    comp.name.padj <- sprintf("%s_padj", comp.name)
    comp.name.fc <- sprintf("%s_log2FoldChange", comp.name)
    expt <- strsplit(split = "_vs_", comp.name)[[1]][1]
    ctrl <- strsplit(split = "_vs_", comp.name)[[1]][2]
    comp.name.plot <- paste(" (",expt,"/",ctrl, ")", sep = "")

    ##Clean up df.diff to remove cases with padj is NA and make 0 as -10^-100
    df.diff <- df.diff[which(!is.na(df.diff[[comp.name.padj]])), ]
    df.diff[[comp.name.padj]] <- ifelse(df.diff[[comp.name.padj]] == 0, 10^-100, df.diff[[comp.name.padj]])

    #Create -log10(padj) column and limit it on thresholds
    df.diff[["log10padj"]] <- -log10(df.diff[[comp.name.padj]])
    if(is.null(pval.plot.cutoff)){
        pval.plot.cutoff <- max(df.diff$log10padj)
    }
    df.diff$log10padj <- ifelse(df.diff$log10padj >= pval.plot.cutoff, pval.plot.cutoff, df.diff$log10padj)
    df.diff$log10padj <- as.numeric(df.diff$log10padj)

    ##Limit foldchange column on thresholds
    if(is.null(FC.plot.cutoff)){
        FC.plot.cutoff <- ceiling(max(abs(range(df.diff[[comp.name.fc]]))))
    }
    df.diff[[comp.name.fc]] <- ifelse(df.diff[[comp.name.fc]] < -FC.plot.cutoff, -FC.plot.cutoff, df.diff[[comp.name.fc]])
    df.diff[[comp.name.fc]] <- ifelse(df.diff[[comp.name.fc]] > FC.plot.cutoff, FC.plot.cutoff, df.diff[[comp.name.fc]])
    df.diff[[comp.name.fc]] <- as.numeric(df.diff[[comp.name.fc]])

    df.diff <- df.diff[ ,c("gene_id", "symbol", comp.name.fc, "log10padj")]
    #Create -log10(padj) column and limit it on thresholds
    if(is.null(go.id)){
    	go.gene.ids <- df.diff$gene_id
    }

    df.diff.go.sig <- subset(df.diff, gene_id %in% go.gene.ids & log10padj >= log10adj.cutoff)
    df.diff.go.insig <- subset(df.diff, gene_id %in% go.gene.ids & log10padj < log10adj.cutoff)
    if(label){
    	df.diff.go.label <- df.diff.go.sig[which(abs(df.diff.go.sig[[comp.name.fc]]) >= log2(label.fc.cutoff)), ]
    	if(label.log10adj.cutoff != log10adj.cutoff){
        df.diff.go.label <- df.diff.go.sig[which((abs(df.diff.go.sig[[comp.name.fc]]) >= log2(label.fc.cutoff)) & df.diff.go.sig$log10padj >= label.log10adj.cutoff), ]
        df.diff.go.unlabel1 <- df.diff.go.sig[which((abs(df.diff.go.sig[[comp.name.fc]]) < log2(label.fc.cutoff)) & df.diff.go.sig$log10padj >= label.log10adj.cutoff), ]
        df.diff.go.unlabel2 <- df.diff.go.sig[which((abs(df.diff.go.sig[[comp.name.fc]]) < log2(label.fc.cutoff)) & df.diff.go.sig$log10padj < label.log10adj.cutoff), ]
        df.diff.go.unlabel3 <- df.diff.go.sig[which((abs(df.diff.go.sig[[comp.name.fc]]) >= log2(label.fc.cutoff)) & df.diff.go.sig$log10padj < label.log10adj.cutoff), ]
        df.diff.go.unlabel <- rbind(df.diff.go.unlabel1, df.diff.go.unlabel2, df.diff.go.unlabel3)
      } else {
        df.diff.go.unlabel <- df.diff.go.sig[which(abs(df.diff.go.sig[[comp.name.fc]]) < log2(label.fc.cutoff)), ]
      }
    }
    df.diff.not.go <- subset(df.diff, !(gene_id %in% go.gene.ids))

    ##generate p with plotting characteristics
    p <- ggplot()
    if(nrow(df.diff.not.go) > 0){
        p <- p + geom_point(data = df.diff.not.go, aes_string(x = comp.name.fc, y = "log10padj", alpha = 0.6), size = point.size, colour = "grey")
    }
    if(nrow(df.diff.go.insig) > 0){
        p <- p + geom_point(data = df.diff.go.insig, aes_string(x = comp.name.fc, y = "log10padj"), size = point.size, colour = "#03467c", alpha = 0.4)
    }
    if(label){
    	if(nrow(df.diff.go.unlabel) > 0){
    	    p <- p + geom_point(data = df.diff.go.unlabel, aes_string(x = comp.name.fc, y = "log10padj"), size = point.size, colour = "#03467c")
    	}
    	if(nrow(df.diff.go.label) > 0){
    	    p <- p + geom_point(data = df.diff.go.label, aes_string(x = comp.name.fc, y = "log10padj"), size = point.size, colour = "#8e0204")
    	}    	
    } else {
    	df.diff.go.sig.fc <- df.diff.go.sig[which(abs(df.diff.go.sig[[comp.name.fc]]) >= log2(fc.cutoff)), ]
    	df.diff.go.sig.nofc <- df.diff.go.sig[which(abs(df.diff.go.sig[[comp.name.fc]]) < log2(fc.cutoff)), ]
    	p <- p + geom_point(data = df.diff.go.sig.nofc, aes_string(x = comp.name.fc, y = "log10padj"), size = point.size, colour = "#03467c")
    	p <- p + geom_point(data = df.diff.go.sig.fc, aes_string(x = comp.name.fc, y = "log10padj"), size = point.size, colour = "#8e0204")
    }


    p <- p + geom_hline(yintercept=log10adj.cutoff, linetype = "dashed")
    p <- p + xlim(c(-FC.plot.cutoff, FC.plot.cutoff))
    p <- p + ylim(c(0, pval.plot.cutoff))
    p <- p + ggtitle(sprintf("%s", go.id))
    p <- p + ylab(expression(paste("-", log[10], " (padjust)", sep = "")))
    p <- p + xlab(substitute(paste(log[2],nn, sep = ""), list(nn = comp.name.plot)))
    if(label && (nrow(df.diff.go.label) > 0)){
        p <- p  + geom_text(data = df.diff.go.label, mapping = aes_string(x = comp.name.fc, y = "log10padj", label = "symbol"),
        hjust = 0.5, vjust = -0.75, fontface = 1, size = label.size)
    }
    p <- p + theme_bw()
    p <- p + theme(plot.title = element_text(size = rel(0.8)),
          axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16), 
          axis.title.x = element_text(size = 18, vjust = -0.5),
          axis.title.y = element_text(size = 18, vjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1, 
          legend.key = element_blank(), legend.text = element_text(size = 14),
          legend.title = element_blank(), legend.position = "none")

    print(p)

}

plotMA_local <- function(df.diff, pvalue=0.05, ylim=NULL, y.qtile=0.98, 
	comp.name=NULL,
    linecol="#ff000080",
    xlab = "mean of normalized counts", #ylab = expression(log[2]~fold~change),
    log.x=TRUE, cex=1, legend.cex = 1, legend.bg = "#FFFFFF88", ...) {

  	comp.name.padj <- sprintf("%s_padj", comp.name)
  	comp.name.fc <- sprintf("%s_log2FoldChange", comp.name)
  	comp.name.bm <- sprintf("%s_baseMean", comp.name)
  	expt <- strsplit(split = "_vs_", comp.name)[[1]][1]
  	ctrl <- strsplit(split = "_vs_", comp.name)[[1]][2]
  	comp.name.plot <- paste(" (",expt,"/",ctrl, ")", sep = "")
	
  	y.axis <- comp.name.fc
  	x.axis <- comp.name.bm
  	ylab <- y.axis
	
  	df.diff <- df.diff[!is.na(df.diff[[comp.name.padj]]) & df.diff[[x.axis]] != 0, ]

	col <- ifelse(df.diff[[comp.name.padj]] <= pvalue & df.diff[[comp.name.fc]] < 0, "#03467c", "#a3a1a133")
	col <- ifelse(df.diff[[comp.name.padj]] <= pvalue & df.diff[[comp.name.fc]] > 0, "#8e0204", col)
	Nhigh <- length(which(col == "#8e0204"))
	Nlow <- length(which(col == "#03467c"))
  	py <- df.diff[[y.axis]]
  	xx <- df.diff[[x.axis]]
  	if (log.x) {
  	  xx <- log2(xx)
  	}
	
  	if (is.null(ylim)) {
  	  ylim <- c(-1,1) * quantile(abs(py[is.finite(py)]), probs=y.qtile) * 1.1
  	}
	
  	pch <- ifelse(py < ylim[1], 6, ifelse(py > ylim[2], 2, 16))
  	py[py <= ylim[1]] <- jitter(rep(ylim[1], sum(py <= ylim[1])))
  	py[py >= ylim[2]] <- jitter(rep(ylim[2], sum(py >= ylim[2])))
	
  	plot(xx, py, cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, pch=pch, 
  	  font =2, font.lab = 2,...)
  	abline(h=0, lwd=4, col=linecol)

    legend("bottomright", text.col = c("black", "#8e0204", "#03467c"), 
    	legend = c(sprintf("Total (%d)", nrow(df.diff)),
        paste("Higher in", expt, sprintf("(%d)", Nhigh)),
        paste("Lower in", expt, sprintf("(%d)", Nlow))),
        cex = legend.cex, bg = legend.bg)	
  	invisible(cbind(data.frame(x=xx, y=py), df.diff))

}

findGeneIdsForGOIds <- function(go.ids){
    xx <- as.list(org.Hs.egGO2ALLEGS)
    gene.ids <- unique(unlist(mclapply(go.ids, function(go.id) xx[[go.id]], mc.cores= 8)))
}

findAllGOOffspring <- function(id){
	xx <- as.list(GOBPOFFSPRING)
	xx <- xx[!is.na(xx)]
	xx[[id]]
}

########################################################################
##This function calculats the rpkm/tpm/tpbm(or rpbm) of the transcripts
##Default is rpkm
##If each of the input SE objects contain the library size then it can
##use that for calculating the rpkm/tpm/tpbm
##For rpkm and tpbm, it finds the length of the transcripts using GRanges
##annotation object then it sums over the exons to give the total number 
##of the reads mapping to the transcript
##Returns  a data frame with rpkm values
########################################################################
RPKM <- function(se, tpbm = F, tpm = F, from.se = TRUE){

    m.gene.count <- assays(se)$count
    multiply.by <- 10^9

    if(from.se){
      vec.library.size <- colData(se)$library.size
      } else {
        vec.library.size <- colSums(m.gene.count)
      }
    
    ##Get the length of the transcripts based on the type ofrow data
    if(inherits(rowRanges(se), "GRangesList")){
      vec.transcript.length <- sum(width(rowRanges(se)))
      } else if(inherits(rowRanges(se), "GRanges")) {
        vec.transcript.length <- width(rowRanges(se))
      }

    if(tpbm || tpm){
      multiply.by <- 10^6
    }

    ##Do not proceede if there is problem with the number of rows
    stopifnot(length(vec.transcript.length) == nrow(m.gene.count))
    if(!tpm){
      m.rpkm <- as.matrix(apply(m.gene.count, 2, function(x){
        x/vec.transcript.length
        }))
      if(nrow(m.gene.count) == 1) m.rpkm <- t(m.rpkm)
    } else {
      m.rpkm <- m.gene.count
    }
    m.rpkm <- t(apply(m.rpkm, 1, function(x){
      multiply.by * x/vec.library.size
      }))
  
    return(m.rpkm)
}


doMyPCA_Heatmap <- function(se, assay.name = NULL, doLog2 = TRUE, main.sub = NULL, 
  doMAD = FALSE, pFeatures = 0.5, scale.type = "none", col.pal = NULL, 
  given.levels = NULL, given.order.cols = NULL, condition.clmn = "condition"){

  if(doLog2){
    m.tpm <- (log2(assays(se)[[assay.name]] + 1))
    } else {
    m.tpm <- (assays(se)[[assay.name]])
  }

  check.na <- apply(m.tpm, 1, function(x) any(is.na(x)))
  m.tpm <- m.tpm[!check.na, ]

  all.var <- apply(t(m.tpm), 2, var)
  if(length(which(all.var == 0)) > 0){
    m.tpm <- m.tpm[-which(all.var == 0),]
  }

  if(doMAD){
    mads <- apply(m.tpm, 1, mad)
    l <- floor(length(mads)*pFeatures)
    m.tpm <- m.tpm[order(mads, decreasing = T)[1:l], ]
  }

  main.pca <- sprintf("PCA on %d (%s)", nrow(m.tpm), assay.name)
  main.mds <- sprintf("MDS on %d (%s)", nrow(m.tpm), assay.name)
  main.heatmap <- sprintf("Heatmap on %d (%s)", nrow(m.tpm), assay.name)
  
  if(!is.null(main.sub)){
    main.pca <- sprintf("%s: %s", main.pca, main.sub)
    main.mds <- sprintf("%s: %s", main.mds, main.sub)
    main.heatmap <- sprintf("%s: %s", main.heatmap, main.sub)
  }

  pca <- prcomp(t(m.tpm), retx = T, scale = T, center = T)
  m.summary <- summary(pca)$importance
  print(m.summary)
  pc1.var <- m.summary[2, 1] * 100
  pc2.var <- m.summary[2, 2] * 100
  df <- as.data.frame(pca$x)
  df$samples <- rownames(df)
  #plot(df$PC1, df$PC2, main = main.pca,
  #    xlab=sprintf("PC1 (%.2f%% of variance)", pc1.var), 
  #    ylab=sprintf("PC2 (%.2f%% of variance)", pc2.var), 
  #    pch=18, col="blue", font = 2, font.lab = 2, cex.main = 0.8)
  #text(df$PC1, df$PC2, row.names(df), cex=0.8, pos=4, col="red", font = 2)

  if(condition.clmn %in% names(colData(se))){
    df$condition <- colData(se[ ,rownames(df)])[[condition.clmn]]
  }

  print(ggplot(df, aes(x = PC1, y = PC2, label = samples, colour = condition)) + 
  geom_point(size = 2) + 
  xlab(sprintf("PC1 (%.2f%% of variance)", pc1.var)) +
  ylab(sprintf("PC2 (%.2f%% of variance)", pc2.var)) +
  ggtitle(main.pca) +
  geom_text(hjust= 0.5, vjust= -0.75, size = 2) +
  theme_bw() + theme(plot.title = element_text(size = rel(0.8)),
  axis.text.x = element_text(size = 16, angle = 60, hjust = 1),
  axis.text.y = element_text(size = 16), 
  axis.title.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  aspect.ratio = 1))

########################################################
#d <- dist(t(m.tpm))
#fit <- cmdscale(d,eig=TRUE, k=2)

## plot solution
#x <- fit$points[,1]
#y <- fit$points[,2]
#xval <- max(abs(max(x)), abs(min(x)))
#yval <- max(abs(max(y)), abs(min(y)))
#df <- data.frame(x, y)
#df$samples <- rownames(fit$points)
#if(condition.clmn %in% names(colData(se))){
#  df$condition <- colData(se[ ,rownames(df)])[[condition.clmn]]
#}
#lim.val <- max(xval, yval)
#plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main = main.mds, 
#  type="n", font = 2, font.lab = 2, 
#  xlim = c(-lim.val, lim.val), ylim = c(-lim.val, lim.val), cex.main = 0.8)
#text(x, y, labels = rownames(fit$points), cex=.8, font = 2)

#print(ggplot(df, aes(x, y, label = samples, colour = condition)) + 
#geom_point(size = 2) + 
#xlab("Coordinate 1") +
#ylab("Coordinate 2") +
#xlim(c(-lim.val, lim.val)) +
#ylim(c(-lim.val, lim.val)) +
#ggtitle(main.mds) +
#geom_text(hjust= 0.5, vjust= -0.75, size = 2) +
#theme_bw() + theme(plot.title = element_text(size = rel(0.8)),
#axis.text.x = element_text(size = 16, angle = 60, hjust = 1),
#axis.text.y = element_text(size = 16), 
#axis.title.x = element_text(size = 18),
#axis.title.y = element_text(size = 18),
#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#aspect.ratio = 1))

###########################################
  m.tpm <- (assays(se)[[assay.name]])

  check.na <- apply(m.tpm, 1, function(x) any(is.na(x)))
  m.tpm <- m.tpm[!check.na, ]

  m.tpm <- normalizeMatrix(m.tpm, doLog2)

  all.var <- apply(t(m.tpm), 2, var)
  if(length(which(all.var == 0)) > 0){
    m.tpm <- m.tpm[-which(all.var == 0),]
  }

  if(doMAD){
    mads <- apply(m.tpm, 1, mad)
    l <- floor(length(mads)*pFeatures)
    m.tpm <- m.tpm[order(mads, decreasing = T)[1:l], ]
  }


  cluster.order <- hclust(dist(m.tpm))$order
  m.tpm <- m.tpm[cluster.order, ]

  if(is.null(col.pal)){
      col.pal <- c("#2F568B", "#ffffff", "#A41E21")
  }

  if( !is.null(given.order.cols) & (condition.clmn %in% names(colData(se)))){
      annotation <- data.frame(sample = factor(colData(se[ ,given.order.cols])[[condition.clmn]]), check.names = F)
      if(!is.null(given.levels)){
        annotation$sample <- factor(annotation$sample, levels = given.levels)
      }
      aheatmap(m.tpm[ , given.order.cols] , Rowv = NA, Colv = NA, 
        main = main.heatmap, legend = T,
        annCol = annotation,
        cexRow = 0.1, col = col.pal,
        scale = scale.type, cex = 0.8)

    } else if(is.null(given.order.cols) & ("condition" %in% names(colData(se)))){
      annotation <- data.frame(sample = factor(colData(se)$condition), check.names = F)
      aheatmap(m.tpm , Rowv = NA, Colv = TRUE, 
        main = main.heatmap, legend = T,
        annCol = annotation, col = col.pal,
        cexRow = 0.5, 
        scale = scale.type, cex = 0.8)
    } else if(is.null(given.order.cols) & !("condition" %in% names(colData(se)))){
      aheatmap(m.tpm , Rowv = NA, Colv = TRUE, 
        main = main.heatmap, legend = T,
        cexRow = 0.1, col = col.pal,
        scale = scale.type, cex = 0.8)
    }
}

##Quantile normalize the features
normalizeMatrix <- function(m, doLog2 = T){

  if(doLog2){
    ##Log2 of the values
    cu <- log2(m + 1)
  } else {
    cu <- m
  }

  ##Quantile normalize the values
  cu.qq <- normalize.quantiles(cu)
  colnames(cu.qq) <- colnames(m)
  rownames(cu.qq) <- rownames(m)

  return(cu.qq)

}

###################
plotReplicates <- function(se.atlas, main = ""){

  df.rpkm <- as.data.frame(log2(assays(se.atlas)$rpkm + 1))
  min.axis.val <- floor(min(df.rpkm))
  max.axis.val <- ceiling(max(df.rpkm))

  samples <- colnames(df.rpkm)
  m.combinations <- combn(samples, m = 2)
  n.combinations <- ncol(m.combinations)
  for(idx in seq(1, n.combinations)){
    xsample.name <- m.combinations[1 ,idx]
    ysample.name <- m.combinations[2 ,idx]
    cor.txt <- sprintf("%d expressed genes, r: %.2f", nrow(df.rpkm), 
      cor(df.rpkm[[xsample.name]], df.rpkm[[ysample.name]]))
    print(ggplot(df.rpkm, aes_string(x = xsample.name, y = ysample.name)) +
      geom_point(size = 2, colour = "#898989") +
      geom_abline(intercept=0, slope = 1) +  
      xlab(bquote(log[2](.(xsample.name)~FPKM~+1))) +
      ylab(bquote(log[2](.(ysample.name)~FPKM~+1))) +
      xlim(c(min.axis.val, max.axis.val)) +
      ylim(c(min.axis.val, max.axis.val)) +
      geom_text(x = (min.axis.val + 9), y = max.axis.val , label = cor.txt, vjust=1, hjust=1, 
        fontface = 1, size = 6) +
      theme_bw() + 
      ggtitle(main) +
      theme(plot.title = element_text(size = rel(0.8)),
          axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16), 
          axis.title.x = element_text(size = 18, vjust = -0.5),
          axis.title.y = element_text(size = 18, vjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1, 
          legend.key = element_blank(), legend.text = element_text(size = 14),
          legend.title = element_blank()))

  }
    
}
