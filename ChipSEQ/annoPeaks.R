# Annotate ChIPseq peaks from Kim Aldinger

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(ggplot2)
#chipDir <- "/.mounts/labs/pailab/private/xsun/output/ncMutMB/20240314/cre/Aldinger-FetalCB"

chipDir <- "/home/rstudio/isilon/private/xsun/output/ncMutMB/20240314/cre/Aldinger-FetalCB/raw"
geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v42.basic.annotation.gtf"
outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_FetalCBL_ChipSeq/output"

dt <- format(Sys.time(), "%Y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}
logFile <- sprintf("%s/annoPeaks.log", outDir)
sink(logFile,  split=TRUE)
tryCatch({
cat("Importing ChIP-seq peaks...\n")
files <- list(
    h3k27ac_1 = sprintf("%s/1_27907-102M_H3K27Ac_hg38.bed", chipDir),
    h3k27ac_2 = sprintf("%s/2_27556-132M_H3K27Ac_hg38.bed", chipDir),
    h3k4me3_1 = sprintf("%s/3_27907-102M_H3K4me3_hg38.bed", chipDir),
    h3k4me3_2 = sprintf("%s/4_27556-132M_H3K4me3_hg38.bed", chipDir),
    pol2_1 = sprintf("%s/5_27907-102M_TotalPol2_hg38.bed", chipDir),
    pol2_2 = sprintf("%s/6_27556-132M_TotalPol2_hg38.bed", chipDir)
)

peaks <- list()
for (nm in names(files)){
    cat("Reading peaks for file:", nm, "\n")
    peaks[[nm]] <- readPeakFile(files[[nm]])
}

# count peaks and print
cat("Counting peaks...\n")
for (nm in names(peaks)) {
    cat(sprintf("Number of peaks for %s: %d\n", nm, length(peaks[[nm]])))
}

#### coverage plot
###cat("Plotting coverage for peaks...\n")
###for (nm in names(peaks)) {
###    cat("Plotting coverage for:", nm, "\n")
###    pdfFile <- sprintf("%s/%s_coverage.pdf", outDir, nm)
###    p <- covplot(peaks[[nm]])
###    #p <- p + theme_minimal(base_size = 12)
###    ggsave(p, file=pdfFile)
###}


#### importing gene definition
###  ### load GENCODE annotation ###
###message("Reading in GENCODE annotation from: %s", geneFile)
###genes <- rtracklayer::readGFF(geneFile)
###genes <- subset(genes, gene_type %in% "protein_coding" & type == "gene")
###  
###geneGR <- GRanges(genes$seqid, 
###                    IRanges(genes$start, genes$end),
###                    name=genes$gene_name, strand=genes$strand
###) 
###promoter <- promoters(geneGR, upstream=2000, downstream=2000)

# Annotate peaks with ChIPseeker
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
cat("Annotating peaks with ChIPseeker...\n")
tagMatrix <- list()
for (nm in names(peaks)) {
    cat("Annotating peaks for:", nm, "\n")
   tagMatrix[[nm]] <- getTagMatrix(peaks[[nm]], windows=promoter)
}

# Plotting tag matrix heatmap
for (nm in names(tagMatrix)) {
    cat("Plotting tag matrix heatmap for:", nm, "\n")
    pdfFile <- sprintf("%s/%s_avgProfileTSS.pdf", outDir, nm)
    p <- plotAvgProf(tagMatrix[[nm]], 
        xlab="Distance to TSS (bp)", 
        ylab="Read Count", 
        xlim=c(-3000, 3000))
    p <- p + ggtitle(sprintf("Tag Matrix Heatmap for %s", nm))
    p <- p + theme_minimal(base_size = 20)
    ggsave(p, file=pdfFile)

    # plot heatmap
    cat("plotting heatmap\n")
    pngFile <- sprintf("%s/%s_tagMatrix_heatmap.png", outDir, nm)
    p <- tagHeatmap(tagMatrix[[nm]])
    p <- p + ggtitle(sprintf("Tag Matrix Heatmap for %s", nm))
    p <- p + theme_minimal(base_size = 20)
    ggsave(p, file=pngFile, width=8, height=24)
    cat("..done\n")
}

peakAnno <- list()
# Annotate peaks with genomic features
for (nm in names(peaks)) {
    cat("Annotating peaks for:", nm, "\n")
    peakAnno[[nm]] <- annotatePeak(peaks[[nm]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
}
# Plotting annotation pie chart
for (nm in names(peaks)) {
    cat("Annotating peaks for:", nm, "\n")
    cat("Plotting annotation pie chart for:", nm, "\n")
    pdfFile <- sprintf("%s/%s_annotation_pie_chart.pdf", outDir, nm)
    pdf(file=pdfFile, width=10, height=8)
    plotAnnoPie(peakAnno[[nm]], 
        main=sprintf("Annotation Pie Chart for %s", nm))
    dev.off()
}

# joint average profile
cat("Plotting joint average profile...\n")
p <- plotAvgProf(tagMatrix, 
    xlab="Distance to TSS (bp)", 
    ylab="Read Count", 
    xlim=c(-3000, 3000))
p <- p + ggtitle("Joint Average Profile for All Peaks")
p <- p + theme_minimal(base_size = 20)
pdfFile <- sprintf("%s/joint_avgProfileTSS.pdf", outDir)
ggsave(p, file=pdfFile, width=10,height=8)

}, error=function(e) {
    cat("An error occurred: ", e$message, "\n")
}, finally = {
    sink()  # Close the log file
})