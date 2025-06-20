## Wrapper for functions for doing gene ontology enrichment tests using
## GOstats
library(AnnotationDbi)
library(annotate)
library(GOstats)
library(xtable)

##
## Help with GO in R:
##   http://www.economia.unimi.it/projects/marray/2007/material/day4/Lecture7.pdf
##
## For humans, use annotation="org.Hs.eg.db" and use entrez ids
## For yeast, use annotation='org.Sc.Sgd', and use ORF ids

##' Runs a GO report for all ontologies
##'
##' The GO tests are run and one html file is written that has the results
##' for all ontologies tested
##'
##' @param genes The entrez ids of your "special" genes
##' @param universe The entrez ids of the "universe"
##' @return The list of results from GOstats::hyperGTest
do.goReport <- function(genes, universe, ontologies=c('BP', 'MF', 'CC'),
                        conditional=TRUE, p.value=0.05, testDirection='over',
                        annotation='org.Hs.eg.db', path=".",
                        report.name="GO-report.html", verbose=FALSE) {
  genes <- unique(as.character(genes))
  universe <- union(genes, as.character(universe))

  gos <- lapply(ontologies, function(ontology) {
    if (verbose) {
      cat("...", ontology, "\n")
    }
    do.gostats(genes, universe, conditional=conditional, p.value=p.value,
               ontology=ontology, annotation=annotation,
               testDirection=testDirection)
  })
  names(gos) <- ontologies

  if (length(grep("\\.html?$", report.name)) == 0L) {
    report.name <- paste(report.name, 'html', sep='.')
  }

  html.out <- file.path(path, report.name)
  if (verbose) {
    cat("Creating HTML report to:", html.out, "\n")
  }

  for (ontology in ontologies) {
    append <- ontology != ontologies[1]
    go.htmlReport(gos[[ontology]], append=append, file=html.out)
    return(gos[[ontology]])
  }

  # return(gos)
  invisible(gos)
}

################################################################################
## Functions below work with one report at a time
## -----------------------------------------------------------------------------
##' Performs GO analysis on a list of selected genes, given the universe.
##' @param genes A list of the genes that are "picked". Expected to be entrez
##' id's for everything except yeast -- in which case, use ORFs
##' @param universe A list of genes that could have been picked.
do.gostats <- function(genes, universe, conditional=TRUE, p.value=0.05,
                       ontology=c('BP', 'MF', 'CC'), annotation='org.Hs.eg.db',
                       testDirection='over') {
  if (!require(GOstats)) {
    stop("This function requires GOstats")
  }
  if (missing(universe)) {
    stop("Gimme the universe of genes we picked from")
  }
  ontology <- match.arg(ontology)
  genes <- as.character(genes)
  universe <- as.character(universe)
  params <- new("GOHyperGParams", geneIds=genes,
                universeGeneIds=universe, annotation=annotation,
                ontology=ontology, pvalueCutoff=p.value,
                conditional=conditional, testDirection=testDirection)
  GO <- hyperGTest(params)
  GO
}


##' Get the genes responsible for GO enrichment in each category. This is a less-robust
##' version of the GOstats::probeSetSummary
##'
##' It will return you the names of the genes in each category, as they were passed
##' in to the do.gostats function. We assume that the genes were passed to the
##' do.gostats function, and that these id's were entrez-like (in yeast, the default
##' is to use the ORF.
##'
##' Returns a data.frame for each enriched GO term. Each data.frame has the
##' entrez-like id, as well as its symobl (if found)
##'
##' @param result A GOHyperGResult object
##' @param p.value The value to use as a pvalue cutoff for the GO categorys
##' @param categorySize Not used
##' @param selected Really, don't use it.
go.members <- function(result, p.value=pvalueCutoff(result), categorySize=NULL,
                       selected=geneIds(NULL)) {
  if (!require('annotate')) {
    stop("Annotate library required")
  }
  if (!is(result, "GOHyperGResult")) {
    stop("result must be a GOHyperGResult instance (or subclass)")
  }

  ## Was this necessary for S. cerevisae?
  ## elements2entrez <- getAnnMap("ENTREZID", annotation(result))
  ## elementIds <- ls(elements2entrez)
  summary <- Category:::XXX_getSummaryGeneric_XXX()
  goids <- summary(result, p.value, categorySize)[, 1] ## Enriched go ids

  universe <- geneIdUniverse(result)[goids]
  sig.ids <- geneIds(result) ## The "interesting" genes from original call

  sig.in.go <- lapply(universe, function(ids) {
    ids <- as.character(ids)
    have <- ids %in% sig.ids
    ids[have]
  })

  if (annotation(result) == 'org.Sc.sgd') {
    name.map <- getAnnMap("GENENAME", annotation(result))
  } else {
    name.map <- getAnnMap("SYMBOL", annotation(result))
  }

  df <- lapply(sig.in.go, function(x) {
    # f <- selectMethod("mget", c("character", "AnnDbBimap"));
    # symbols <- f(x, name.map, ifnotfound=NA)
    symbols <- AnnotationDbi::mget(x, name.map, ifnotfound=NA)
    symbols <- sapply(symbols, '[', 1)
    symbols[is.na(symbols)] <- names(symbols)[is.na(symbols)]
    d.f <- data.frame(id=x, symbol=symbols)
    rownames(d.f) <- NULL
    d.f
  })

  df
}

##' An augmented HTML report for the enriched GO terms.
##'
##' This includes a column where the elements found for the enrichment are
##' included.
##'
##' @param result A GOHyperGResult object
##' @param with.elements Whether or not (TRUE/FALSE) to print the genes that
##' were used to determine enrichment for the particular category.
go.htmlReport <- function(result, p.value=pvalueCutoff(result), file="",
                          append=FALSE, label="", digits=3,
                          htmlLinks=TRUE, with.elements=TRUE,
                          link.elements=TRUE) {
  have.xtable <- suppressWarnings({
    require("xtable", quietly=TRUE, warn.conflicts=FALSE)
  })
  if (!have.xtable) {
    stop("xtable required for htmlReport", call.=FALSE)
  }
  if (!is(result, "GOHyperGResult")) {
    stop("result must be a GOHyperGResult instance (or subclass)")
  }
  if (nchar(file) == 0) {
    file <- sprintf("enriched.%s.p%.2f.html",
                    paste(testName(result), collapse=","),
                    p.value)
  }

  AMIGO_URL <- paste("http://www.godatabase.org/cgi-bin/amigo/go.cgi?",
                     "view=details&search_constraint=terms&depth=0&query=%s",
                     sep="")

  ## Figure out the URL to use for element linking
  organism <- unlist(strsplit(annotation(result), '\\.'))[2]
  if (organism == 'Sc') {
    ELEMENT_URL <- "http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=%s"
  } else {
    ELEMENT_URL <- "http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s"
  }

  summary <- Category:::XXX_getSummaryGeneric_XXX()
  df <- do.call(summary, list(result, p.value, htmlLinks=htmlLinks))

  if (nrow(df) == 0) {
    warning("No results to report")
    return(invisible(TRUE))
  }

  if (with.elements) {
    members <- go.members(result, p.value)
    df$members <- lapply(df[,1], function(go.id) {
      elems <- members[[go.id]]$symbol

      if (length(elems) > 0) {
        if (link.elements) {
          urls <- sprintf(ELEMENT_URL, elems)
          elems <- sprintf('<a href="%s">%s</a>', urls, elems)
        }
        elems <- paste(elems, collapse=", ")
      } else {
        warning("Empty symbol list for:", go.id)
        elems <- ""
      }
      elems
    })

  }

  dig <- rep(digits, ncol(df) + 1) ## need +1 for xtable row name
  dig[5:7] <- 0
  xt <- xtable(df, caption=paste(testName(result), collapse="/"), digits=dig)
  print(xt, type="html", file=file, append=append, caption.placement='top',
        sanitize.text.function=identity, include.rownames=FALSE)
  invisible(list(xtable=xt, df=df))
}
