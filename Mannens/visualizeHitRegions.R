# use GViz to examine regions where DMRs overlap with Cicero connections and TSS sites.
rm(list=ls())
library(Gviz)
library(cicero)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Organism.dplyr)
library(org.Hs.eg.db)

inFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/Mannens_Integration/250612/Mannens_DMR_conns_annotated_250612.txt"
dmrFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"

mannensDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024_Cicero/250611"
ciceroFile <- sprintf("%s/Mannens_cicero_conns_allchroms.qs",mannensDir)

geneDef <- "/home/rstudio/isilon/src/gencode/GRCh38/gencode.v42.basic.annotation.gtf"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/GvizHitRegions"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

cat("* Reading DMR connections from file\n")
dmr_conns <- read.delim(inFile, sep="\t", h=T)
cat(sprintf("Read %i DMR connections from file\n", nrow(dmr_conns)))

dmr <- read.delim(dmrFile, sep="\t", h=T)
dmr_GR <- makeGRangesFromDataFrame(dmr, 
    seqnames.field = "chr", start.field = "start", end.field = "end",
    keep.extra.columns = TRUE)

gene_anno <- rtracklayer::import(geneDef)
tss <- start(gene_anno)
tss[which(gene_anno$strand == "-")] <- end(gene_anno)[which(gene_anno$strand == "-")]
gene_anno2 <- GRanges(
    seqnames = seqnames(gene_anno),
    ranges = IRanges(start=tss, end=tss),
    gene = gene_anno$gene_name
)

# Convert peak1 into a GenomicRanges object. 
# Find the peaks that overlap with the TSS sites.
g <- dmr_conns
pk1 <- unlist(strsplit(g$Peak1, "-"))
    chrom <- pk1[seq(1,length(pk1), by=3)]
    starts <- as.numeric(pk1[seq(2,length(pk1), by=3)])
    ends <- as.numeric(pk1[seq(3,length(pk1), by=3)])
tmp <- GRanges(chrom, IRanges(starts,ends),name=g$Peak1)
ol <- findOverlaps(tmp,gene_anno2)
both <- as.data.frame(
    cbind(as.character(g$Peak1[queryHits(ol)]), 
          as.character(gene_anno2$gene[subjectHits(ol)]))
)
colnames(both) <- c("Peak1", "gene.Peak1")
both <- both[!duplicated(both),]
cat(sprintf("Found %i peaks overlapping with TSS sites\n", nrow(both)))

merged <- merge(x=dmr_conns, y=both, by="Peak1",all.x=TRUE)

# now do the same for peak2
pk2 <- unlist(strsplit(g$Peak2, "-"))
    chrom <- c(chrom, as.character(pk2[seq(1,length(pk2), by=3)])) 
    starts <- c(starts, as.numeric(pk2[seq(2,length(pk2), by=3)]))
    ends <- c(ends, as.numeric(pk2[seq(3,length(pk2), by=3)]))
tmp <- GRanges(chrom, IRanges(starts,ends),name=g$Peak2)
ol <- findOverlaps(tmp,gene_anno2)
both <- as.data.frame(
    cbind(as.character(g$Peak2[queryHits(ol)]), 
          as.character(gene_anno2$gene[subjectHits(ol)]))
)
colnames(both) <- c("Peak2", "gene.Peak2")
both <- both[!duplicated(both),]
p2 <- merge(x=merged, y=both, by="Peak2",all.x=TRUE)
cat(sprintf("Found %i peaks overlapping with TSS sites\n", nrow(p2)))

rm(gene_anno, gene_anno2, ol, both, tmp, pk1, pk2, chrom, starts, ends)

dmr_conns <- p2
###cat("Reading Cicero connections from file\n")
##conns <- qs::qread(ciceroFile)
##conns[,2] <- as.character(conns[,2])
##blah <- which(conns$Peak1 %in% dmr_conns$Peak1 & conns$Peak2 %in% dmr_conns$Peak2)
##conns <- conns[blah,]
##
# convert Cicero connections to GRanges
g <- dmr_conns
pk1 <- unlist(strsplit(g$Peak1, "-"))
    chrom <- pk1[seq(1,length(pk1), by=3)]
    starts <- as.numeric(pk1[seq(2,length(pk1), by=3)])
    ends <- as.numeric(pk1[seq(3,length(pk1), by=3)])
pk2 <- unlist(strsplit(g$Peak2, "-"))
    chrom <- c(chrom, as.character(pk2[seq(1,length(pk2), by=3)])) 
    starts <- c(starts, as.numeric(pk2[seq(2,length(pk2), by=3)]))
    ends <- c(ends, as.numeric(pk2[seq(3,length(pk2), by=3)]))

df <- data.frame(
    seqnames = chrom,
    start = starts,
    end = ends
)
df <- df[!duplicated(df),]
gr <- makeGRangesFromDataFrame(df, 
    seqnames.field = "seqnames", start.field = "start", end.field = "end",
    keep.extra.columns = TRUE)
rm(pk1, pk2, starts, ends, chrom, g)

cat("Making tracks for plotting\n")
atrack <- GenomeAxisTrack()
dtrack <- AnnotationTrack(range=dmr_GR, name="DMRs",fill="green",col="green")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

pk <- dmr_conns[,c("Peak1","Peak2","coaccess")]
pk <- pk[!duplicated(pk),]

# plot tracks around a ROI using Gviz
plotGviz <- function(chrom,spos,epos,ttl,outFile){
    itrack <- IdeogramTrack(genome="hg38", chromosome=chrom)
    gtrack <- GeneRegionTrack(txdb,
        genome="hg38",
        chromosome=chrom,start=spos,end=epos,
        showId=TRUE, geneSymbol=TRUE
    )
    symbols <- unlist(mapIds(org.Hs.eg.db, gene(gtrack),"SYMBOL","ENTREZID", multiVals="first"))
    symbol(gtrack) <- symbols[gene(gtrack)]
    displayPars(gtrack) <- list(fontsize.group=20)

    ###gr2 <- subsetByOverlaps(gr, GRanges(chrom, IRanges(spos, epos)))
    ###ctrack <- AnnotationTrack(range=gr2, name="Cicero connections",
    ###    col="black", fill="black")

    p <- plot_connections(pk,chr=chrom, minbp=spos, maxbp=epos, 
        alpha_by_coaccess=TRUE,coaccess_cutoff=0.1,
        include_axis_track= FALSE,
        return_as_list=TRUE)

    winwd <- epos - spos
    wd <- min(28,14 * (winwd/20000) * 0.3)
    ht <- min(10,6 * (winwd/20000) * 0.3)

    ###hlight <- HighlightTrack(trackList=list(p[[2]], dtrack, gtrack),
    ###    range=dmr_GR)

    png(outFile, width=wd, height=ht, units="in",res=300)
    tryCatch({
        cat("Calling plotTracks\n")
        plotTracks(
            trackList=list(itrack, atrack, p[[1]],p[[2]],dtrack, gtrack),
            from=spos, to=epos,
            chromosome=chrom,
            sizes=c(0.1, 0.1, 0.2, 0.1, 0.1, 0.5),
            main=ttl
        )
    }, error=function(e) {
        message(sprintf("Error plotting tracks: %s", e$message))
        browser()
    }, finally={
        dev.off()
    })
}

# get a window around a gene of interest
getWins <- function(nm){
    g <- subset(dmr_conns, gene.Peak1 == nm | gene.Peak2 == nm)
    if (nrow(g) == 0) {
        cat(sprintf("No connections found for gene %s\n", nm))
        return(NULL)
    }
    print(head(g))
    pk1 <- unlist(strsplit(g$Peak1, "-"))
    starts <- as.numeric(pk1[seq(2,length(pk1), by=3)])
    ends <- as.numeric(pk1[seq(3,length(pk1), by=3)])
    pk2 <- unlist(strsplit(g$Peak2, "-"))
    starts <- c(starts, as.numeric(pk2[seq(2,length(pk2), by=3)]))
    ends <- c(ends, as.numeric(pk2[seq(3,length(pk2), by=3)]))
    
    spos <- min(starts)
    epos <- max(ends)
    chrom <- as.character(pk1[1])
    ln <- epos - spos

    cat(sprintf("Gene %s: %s %i-%i (%i kb win)\n", 
        nm, chrom, spos, epos, round((epos - spos) / 1000)))

    ###maxwd <- 200*1000
    ###if (ln > maxwd){
    ###    gr <- GRanges(chrom, IRanges(spos, epos))
    ###    gr <- resize(gr, width=maxwd, fix="center")
    ###    spos <- start(gr)
    ###    epos <- end(gr)
    ###    cat(sprintf("Resized to %i kb window\n", round(maxwd / 1000)))
    ###       cat(sprintf("Gene %s: %s %i-%i (%i kb win)\n", 
    ###    nm, chrom, spos, epos, round((epos - spos) / 1000)))
    ###}
    return(list(chrom=chrom, spos=spos, epos=epos))
}

gn <- "GSE1"
wins <- getWins(gn)
if (!is.null(wins)) {
    ttl <- sprintf("%s: %s - %s : %s (%i kb win)",
        wins$chrom, 
        prettyNum(wins$spos,big.mark=",",sci=FALSE),
        prettyNum(wins$epos,big.mark=",",sci=FALSE), 
        gn,
        round((wins$epos - wins$spos) / 1000
    ))
    plotGviz(wins$chrom, wins$spos, wins$epos,ttl=ttl,
        outFile=sprintf("%s/%s.png", outDir, gn))
} else {
    cat("No valid windows found for plotting.\n")
}

###cat("* Plotting...\n")









