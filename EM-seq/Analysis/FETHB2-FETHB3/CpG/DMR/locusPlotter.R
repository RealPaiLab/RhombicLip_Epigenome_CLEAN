# Uses BSseq and BSmooth to plot methylation for a region
library(DSS)

#' Plot methylation for a given region
#' @param bs_mega BSseq object containing methylation data
#' @param region GRanges object or a string in the format "chr-start-end"
#' @param smooWin (integer) Window size for smoothing the methylation data, if NULL no smoothing is applied
#' @param extend Extend view region by this many bp on either side
#' @param cols (char) Colors to use for sample-level lines
#' @param outFile File to which methylation plot is saved
#' @return NULL, saves the plot to a file
methyPlotter <- function(bs_mega, region,extend=2000,                          
                         cols,geneTrack, dmrTrack = NULL,
                         smooWin=1e6,outFile,
                         geneName="gene",
                         lgd_cols = NULL,
                         collapseTranscripts = FALSE) {
  
  # region is a GRanges object
  # up and down are in bp
  if (is.character(region)) {
    x <- unlist(strsplit(region, "-"))
    region <- GRanges(x[1], 
        IRanges(as.numeric(x[2]), as.numeric(x[3])))
  } else if (class(region) == "data.frame"){
    region <- makeGRangesFromDataFrame(region, 
        keep.extra.columns = TRUE, 
        seqnames.field = "chr", 
        start.field = "start", 
        end.field = "end")
  } else if (class(region) != "GRanges") {
    cat("invalid region object\n")
  }

  # apply smoothing if smooWin is provided
    if (!is.null(smooWin)) {
        # extend the region
        region_smoo <- resize(
                region, width = width(region) + smooWin, fix = "center"
        )
        cat("subsetting the object\n")

        bs <- subsetByOverlaps(bs_mega, region_smoo)
        cat("smoothing\n")
        bs_smoo <- BSmooth(bs)
    }   else {
        bs_smoo <- bs_mega
    }
  
  rg <- resize(region, width = width(region) + (extend*2), fix = "center")
  # plot the methylation

bs_smoo <- subsetByOverlaps(bs_smoo, rg)
cat("In methyPlotter\n")
  #pdf(outFile, width=16,height=16)
  #tryCatch({
  ttl <- sprintf("%s: %s %s-%s (%i kb)",
      geneName, seqnames(rg)[1], 
      prettyNum(start(rg), big.mark=","),
      prettyNum(end(rg), big.mark=","),
      round(width(rg)/1000)
    )
  plotRegion_Gviz(bs_smoo,rg,cols,outFile=outFile, ttl=ttl,
    lgd_cols = lgd_cols, dmrTrack = dmrTrack,
    collapseTranscripts = collapseTranscripts)

 ### },error=function(ex){
 ###   print(ex)
 ### }, finally={
 ###   dev.off()
 ### })
}
  
plotRegion_Gviz <- function(BSseq,gr,cols,ttl="gene",outFile="test.pdf",
  lgd_cols = NULL, dmrTrack = NULL, collapseTranscripts = FALSE) {
  library(Gviz)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  ENCODEels <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/GRCh38-ELS.bed"
  els <- read.delim(ENCODEels, sep="\t",header=FALSE, as.is=T); #rtracklayer::import(ENCODEels)
  els <- makeGRangesFromDataFrame(els, 
  keep.extra.columns = TRUE, 
  seqnames.field = "V1", 
  start.field = "V2", 
  end.field = "V3"
)
  els <- subsetByOverlaps(els, gr)
  spos <- start(gr); epos <- end(gr); chrom <- as.character(seqnames(gr))

  if (!is.null(dmrTrack)){
  dmrTrack <- subsetByOverlaps(dmrTrack, gr)
  }
    # now plot the Gviz data
    cat("* Creating Gviz tracks\n")
    cat("* Getting ideogram\n")
    itrack <- IdeogramTrack(genome="hg38", chromosome=chrom,cex=1.5)
    atrack <- GenomeAxisTrack(cex=1.5)
    cat("* Getting gene track\n")
    gtrack <- GeneRegionTrack(txdb,
        genome="hg38",
        chromosome=chrom,start=spos,end=epos,
        showId=TRUE, geneSymbol=TRUE,
        cex = 2.4, lwd=2,
        collapseTranscripts = collapseTranscripts)

    elsTrack <- AnnotationTrack(
        els, 
        name="ELS", 
        genome="hg38", 
        chromosome=chrom,
        fill="red", 
        col="white",
        stacking="dense"
    )

    if (!is.null(dmrTrack)) {
        dmrTrack <- AnnotationTrack(
            dmrTrack, 
            name="DMRs", 
            genome="hg38", 
            chromosome=chrom,
            fill="dodgerblue", 
            col="dodgerblue",
            stacking="dense"
        )
        dmrwd <- 0.1
    } else {
        dmrTrack <- NULL
        dmrwd <- 0
    }

    symbols <- unlist(mapIds(org.Hs.eg.db, gene(gtrack),"SYMBOL","ENTREZID", multiVals="first"))
    symbol(gtrack) <- symbols[gene(gtrack)]
    displayPars(gtrack) <- list(fontsize.group=20)

    winwd <- epos - spos
    wd <- min(28,14 * (winwd/20000) * 0.3)
    ht <- min(10,6 * (winwd/20000) * 0.3)
    if (ht < 8) ht <- 8
    if (wd < 16) wd <- 16

    if (collapseTranscripts != FALSE) {
      geneHt <- 0.3
    } else {
      geneHt <- 0.8
    }

    pdf(outFile, width=wd, height=ht)# units="in",res=300)
    tryCatch({
        cat("Calling plotTracks\n")
        if (is.null(dmrTrack)) {
            plotTracks(
                trackList=list(itrack, atrack, elsTrack, gtrack),
                from=spos, to=epos,
                chromosome=chrom,
                sizes=c(0.1, 0.1, 0.1, geneHt),
                main=ttl
            )
        } else {
            plotTracks(
                trackList=list(itrack, atrack, elsTrack, dmrTrack, gtrack),
                from=spos, to=epos,
                chromosome=chrom,
                sizes=c(0.1, 0.1, 0.1, dmrwd, geneHt),
                main=ttl
            )
        }
        
        cat("Plotting methylation track\n")
        MethyPlot_basic(BSseq, gr, cols,ttl, lgd_cols=lgd_cols,highlightRegion=dmrTrack)
        cat("done\n")
    }, error=function(e) {
        message(sprintf("Error plotting tracks: %s", e$message))
        browser()
    }, finally={
        dev.off()
    })
}  

MethyPlot_basic <- function(BSseq, gr, cols,ttl="gene", lgd_cols = NULL,
  highlightRegion=NULL) {
    positions <- start(BSseq)
    smoothPs <- getMeth(BSseq, type = "smooth")
    rawPs <- getMeth(BSseq, type = "raw")
    coverage <- getCoverage(BSseq)

    
    
    par(xaxs="i")
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
         bty="n",
         ylim = c(0,1), 
         xlim = c(start(gr), end(gr)), 
         xlab = "", ylab = "Methylation",cex.lab=2)         
    axis(side = 2, at = c(0.2, 0.5, 0.8),cex.axis=2)
    axis(side = 1,cex.axis=2)
    rug(positions)
    title(main = ttl, cex.main = 2)

   if (!is.null(highlightRegion)) {
        rect(start(highlightRegion), 0, end(highlightRegion), 1, 
        col=rgb(30/255, 144/255, 1,0.2), border=NA)
    }

    sapply(1:ncol(BSseq), function(sampIdx) {
        .bsPlotSample(positions, smoothPs[, sampIdx], col = cols[sampIdx],
                     lty = 1, lwd = 1,
                     plotRange = c(start(gr), end(gr)))
    })
    legend("bottomright", legend=names(lgd_cols), col=lgd_cols, 
      lty=1, lwd=2, cex=1.5, bty="n")
}

#' plotting function for the Gviz CustomTrack that plots the methylation
#' data
MethyPlot <- function(GdObject, prepare=FALSE) {
    bs_data <- GdObject@variables
    BSseq <- bs_data$BSseq
    gr <- bs_data$gr
    cols <- bs_data$colours

    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames
    positions <- start(BSseq)
    smoothPs <- getMeth(BSseq, type = "smooth")
    rawPs <- getMeth(BSseq, type = "raw")
    coverage <- getCoverage(BSseq)

    ###plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
    ###     ylim = c(0,1), 
    ###     xlim = c(start(gr), end(gr)), 
    ###     xlab = "", ylab = "Methylation")         
    ###axis(side = 2, at = c(0.2, 0.5, 0.8))
    rug(positions)

    sapply(1:ncol(BSseq), function(sampIdx) {
      print(sampIdx)
        .bsPlotSample(positions, smoothPs[, sampIdx], col = cols[sampIdx],
                     lty = 1, lwd = 1,
                     plotRange = c(start(gr), end(gr)))
    })

    return(GdObject)
}

.bsPlotSample <- function(x, y, col, lty, lwd, plotRange) {
    if(sum(!is.na(y)) <= 1)
        return(NULL)
    xx <- seq(from = plotRange[1], to = plotRange[2], length.out = 500)
    yy <- approxfun(x, y)(xx)
    lines(xx, yy, col = col, lty = lty, lwd = lwd)
}

#' Plot DMRs on chromosomes using Gviz
#' @param gr GRanges object containing DMRs
plotDMRsOnChromosomes <- function(gr,outDir){
  # create one plot per chromosome and combine into a single PDF
  library(Gviz)
  cat("* Creating Gviz tracks\n")
  cat("* Getting ideogram\n")
  for (chrom in unique(seqnames(gr)[1:3])) {
    cat(sprintf("Processing chromosome: %s\n", chrom))
    rg <- gr[seqnames(gr) == chrom]
    if (length(rg) == 0) next

    itrack <- IdeogramTrack(genome="hg38", chromosome=chrom)
    atrack <- GenomeAxisTrack()
    
    dtrack <- AnnotationTrack(
      rg, 
      name="DMRs", 
      genome="hg38", 
      chromosome=chrom,
      fill="blue", 
      col="black",
      stacking="dense"
    )

    pdf(sprintf("%s/DMRs_on_%s.pdf", outDir,chrom), width=10, height=6)
    plotTracks(list(itrack, atrack, dtrack), from=start(rg)[1], to=end(rg)[length(rg)])
    dev.off()
  }

  
  


}