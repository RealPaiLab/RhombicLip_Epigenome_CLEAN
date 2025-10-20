# Link DMRs to genes using Mannens dataset

rm(list=ls())
library(DSS)

source("../../utils_PaiLab.R")

cytoDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/CpG_snpFiltered"

phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv"

CPG_DMR_FILE <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"

Mannens_peakLinkFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024/250714/Mannens_linked_peaks_250714.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/Mannens_Integration"

dt <- format(Sys.Date(), "%y%m%d")
outDir <- file.path(outDir, dt)
if (!file.exists(outDir))
    dir.create(outDir, recursive = FALSE, showWarnings = FALSE)

logFile <- sprintf("%s/log_%s.txt", outDir, dt)
sink(logFile,split=TRUE)

tryCatch({
 
cat("* Reading DMRs\n")
cat("DMR file\n")
print(CPG_DMR_FILE, quote = FALSE)

dmrs <- read.delim(
    CPG_DMR_FILE, 
    header = TRUE, sep="\t",
    stringsAsFactors = FALSE
)
dmrs$name <- paste(dmrs$chr, 
    dmrs$start, 
    dmrs$end, 
    sep="-"
)

cat("* Reading Mannens linked peaks\n")
print(Mannens_peakLinkFile, quote = FALSE)
peaks <- read.delim(Mannens_peakLinkFile,sep="\t", header=TRUE, stringsAsFactors = FALSE)

peakGR <- makeGRangesFromDataFrame(peaks, 
    keep.extra.columns = TRUE, 
    seqnames.field = "seqnames", 
    start.field = "start", 
    end.field = "end"
)

dmrGR <- makeGRangesFromDataFrame(dmrs, 
    keep.extra.columns = TRUE, 
    seqnames.field = "chr", 
    start.field = "start", 
    end.field = "end"
)

# find overlaps between DMRs and Mannens peaks. Create a joint data.frame with the DMRs and the Mannens peaks.
overlaps <- findOverlaps(dmrGR, peakGR, ignore.strand = TRUE)
dol <-  as.data.frame(dmrGR[queryHits(overlaps)])
pol <- as.data.frame(peakGR[subjectHits(overlaps)])[,c("peak","gene")]
colnames(pol)[1:2] <- c("Mannens_peak","Mannens_peakGene")
dmrGenes <- cbind(dol, pol)
dmrGenes <- dmrGenes[!duplicated(dmrGenes),]

# add the DMRs that don't overlap with any Mannens peaks, using NA for the Mannens columns
cat("* Adding DMRs without Mannens peaks\n")
miss <- setdiff(1:length(dmrGR), queryHits(overlaps))
dmrMiss <- as.data.frame(dmrGR[miss])
dmrMiss$Mannens_peak <- NA
dmrMiss$Mannens_peakGene <- NA
dmrGenes <- rbind(dmrGenes, dmrMiss)



cat("* Writing output\n")
outFile <- file.path(outDir, "DMRs_linkedToMannensPeaks.csv")
write.table(dmrGenes, file=outFile, sep="\t", row.names = FALSE, quote = FALSE) 

dmrSimpler <- dmrGenes[,-which(colnames(dmrGenes) == "Mannens_peak")]
dmrSimpler <- dmrSimpler[!duplicated(dmrSimpler),]

cat("\n\n-----------------------------------\n\n")
cat("Annotating DMRs\n")
cat("-----------------------------------\n\n")

# print number of unique DMRs and number of unique genes
cat("* Number of unique DMRs: ", length(unique(dmrSimpler$name)), "\n")
cat("* Number of unique Mannens genes: ", length(unique(dmrSimpler$Mannens_peakGene)), "\n\n")

uqGenes <- unique(dmrSimpler$Mannens_peakGene)

  cat("* Getting neurodev genes\n")
  neurodev_genes <- get_neurodevGenes()
  dmr_neurodev <- intersect(uqGenes, neurodev_genes)
  cat(sprintf("%i DMR-associated genes are known neurodevelopmental regulators\n\n",
    length(dmr_neurodev)  ))
  dmrSimpler$IsNeurodev <- dmrSimpler$Mannens_peakGene %in% neurodev_genes
  cat("\nNeurodev genes:\n")
  cat(sprintf("%s\n\n", paste(sort(dmr_neurodev), collapse=", ")))

  cat("* Getting Group 3 & 4 MB associated genes\n")
  g34_genes <- get_g34genes()
    dmr_g34 <- intersect(uqGenes, g34_genes)
    cat(sprintf("%i DMR-associated genes are known Group 3 & 4 MB associated genes\n\n",
        length(dmr_g34)  ))
    dmrSimpler$IsG34gene <- dmrSimpler$Mannens_peakGene %in% g34_genes  
      cat("G34 MB genes:\n")
  cat(sprintf("%s\n\n", paste(sort(dmr_g34), collapse=", ")))

  cat("* Getting human-accelerated regions (HARs)\n") 
  har <- getHARs()
  ol <- findOverlaps(dmrGR, har, ignore.strand = TRUE)  
  dmrSimpler$IsHAR <- FALSE
    if (length(ol) > 0) {
        olDMR <- dmrGR$name[queryHits(ol)]
        dmrSimpler$IsHAR[which(dmrSimpler$name %in% olDMR)] <- TRUE
    }   
   cat(sprintf("%i DMRs overlap with human-accelerated regions (HARs)\n\n",
        sum(dmrSimpler$IsHAR)  ))

    cat("* Getting fetal CBL peaks\n")
    fetalCB_peaks <- getFetalCB_HistonePeaks()
    ol <- findOverlaps(dmrGR, fetalCB_peaks$H3K27ac, ignore.strand = TRUE)
    dmrSimpler$OLFetalCBH3K27ac <- FALSE
    if (length(ol) > 0) {
        olDMR <- dmrGR$name[queryHits(ol)]
        dmrSimpler$OLFetalCBH3K27ac[which(dmrSimpler$name %in% olDMR)] <- TRUE
    }
    cat(sprintf("%i DMRs overlap with fetal CBL peaks\n\n",
        sum(dmrSimpler$OLFetalCBH3K27ac)  ))

cat("Overlap with Northcott 2017 amplifications and deletions\n")
sv <- suppressMessages(getNorthcott2017_AmpsDels())
for (x in c("amps","dels")){
    for (y in names(sv[[x]])) {
        cat(sprintf("%s: %s\n", x, y))
        gr <- sv[[x]][[y]]
        ol <- findOverlaps(dmrGR, gr, ignore.strand = TRUE)
        cur <- rep(FALSE, nrow(dmrSimpler))
        if (length(ol) > 0) {
            olDMR <- dmrGR$name[queryHits(ol)]
            cur[which(dmrSimpler$name %in% olDMR)] <- TRUE
        }
        cat(sprintf("%i DMRs overlap with %s %s\n\n",
            sum(cur), x, y))
        dmrSimpler[,sprintf("OLNorthcott2017_%s_%s", x, y)] <- cur
    }
}

sv <- suppressMessages(getNorthcott2012_AmpsDels())
for (x in c("amps","dels")){
    for (y in names(sv[[x]])) {
        cat(sprintf("%s: %s\n", x, y))
        gr <- sv[[x]][[y]]
        ol <- findOverlaps(dmrGR, gr, ignore.strand = TRUE)
        cur <- rep(FALSE, nrow(dmrSimpler))
        if (length(ol) > 0) {
            olDMR <- dmrGR$name[queryHits(ol)]
            cur[which(dmrSimpler$name %in% olDMR)] <- TRUE
        }
        cat(sprintf("%i DMRs overlap with %s %s\n\n",
            sum(cur), x, y))
        dmrSimpler[,sprintf("OLNorthcott2012_%s_%s", x, y)] <- cur
    }
}

idx <- grep("SHH", colnames(dmrSimpler))
dmrSimpler <- dmrSimpler[,-idx]

# write results to file
outFile <- sprintf("%s/DMRs_linkedToMannensPeaks_annotated_%s.csv",outDir, dt)
write.table(dmrSimpler, file=outFile, sep="\t", 
    row.names = FALSE, quote = FALSE) 

cat("\n\n------------------------------------\n")
cat("Final statistics\n")
cat("------------------------------------\n")
cat(sprintf("Num RL-VZ / RL-SVZ DMRs: %i\n", nrow(dmrs)))
x <- length(unique(dmrSimpler$name[which(!is.na(dmrSimpler$Mannens_peakGene))]))
cat("* Number of unique DMRs linked to Mannens peaks: ", x, "\n")
cat("* Number of unique Mannens genes linked to DMRs: ", length(unique(dmrGenes$Mannens_peakGene)), "\n")
cat("* Number of unique DMRs linked to neurodevelopmental genes: ", 
    length(unique(dmrSimpler$name[which(dmrSimpler[,"IsNeurodev"])])), "\n")
cat(sprintf("* Number of neurodev genes overlapping DMRs = %s\n", 
    length(unique(dmrSimpler$Mannens_peakGene[which(dmrSimpler[,"IsNeurodev"])]))))        
cat("* Number of unique DMRs linked to Group 3 & 4 MB genes: ", 
    length(unique(dmrSimpler$name[which(dmrSimpler[,"IsG34gene"])])), "\n")
cat(sprintf("* Number of Group 3 & 4 MB genes overlapping DMRs = %s\n", 
    length(unique(dmrSimpler$Mannens_peakGene[which(dmrSimpler[,"IsG34gene"])]))))
for (colNum in 17:ncol(dmrSimpler)) {
    cat(sprintf("%s : %s\n", 
        colnames(dmrSimpler)[colNum], 
        length(unique(dmrSimpler$name[which(dmrSimpler[,colNum])])) ))
}

browser()
# Plot DMRs associated with selected genes
cat("Reading methylation data to plot\n")
cat("* Reading BSseq object\n")
t0 <- Sys.time()
  bsObj <- readBS(inDir = cytoDir, 
                  grepPattern = "\\.txt\\.gz$")
print(Sys.time() - t0)

pheno <- read.delim(phenoFile, 
        header = TRUE, 
        stringsAsFactors = FALSE
)
dpos <- regexpr("-",pheno$library_ID)
pheno$sample <- substr(pheno$library_ID, 1, dpos-1)
# define colours for plots
if (all.equal(pheno$sample, sampleNames(bsObj$bs))!=TRUE) {
  stop("Sample names in pheno and bsObj do not match!")
}
cols <- c("#7B3294","#008837")[factor(pheno$ROI,
  levels=c("VZ","SVZ"))]

browser()
inFile <- get_configs("GENCODE_GENE_FILE") 
genes <- rtracklayer::readGFF(inFile)
source("locusPlotter.R")

dmrGR <- makeGRangesFromDataFrame(dmrs, 
    keep.extra.columns = TRUE, 
    seqnames.field = "chr", 
    start.field = "start", 
    end.field = "end"
)

geneList <- "GSE1" #c("EOMES","OTX2","SOX2","WLS")
for (g in geneList){
    print(g)
    curg <- subset(genes, gene_name == g & type == "gene")
    curg$exon_number <- 1
    curg$isoforms <- 1

    methyPlotter(
        bsObj$bs,
        region = getPlottingRange(dmrSimpler, g),        
        geneTrack = curg,
        cols = cols,
        lgd_cols = c("RL-VZ" = "#7B3294", "RL-SVZ" = "#008837"),
        geneName = g,
        outFile = sprintf("%s/%s_allLinkedDMRs.pdf", outDir,g),
        dmrTrack = dmrGR,
        collapseTranscripts = "meta" # see Gviz documentation for GeneRegionTrack
    )
}

regionSet <- list(
    EOMES = GRanges("chr3", IRanges(27716000, 27728599)),
    OTX2 = GRanges("chr14", IRanges(56728210, 56900000)),
    WLS = GRanges("chr1", IRanges(68140000,68210000))
)
for (g in names(regionSet)) {
    print(g)
    methyPlotter(
        bsObj$bs,
        region = regionSet[[g]],        
        geneTrack = genes[genes$gene_name == g & genes$type == "gene",],
        cols = cols,
        lgd_cols = c("RL-VZ" = "#7B3294", "RL-SVZ" = "#008837"),
        geneName = g,
        outFile = sprintf("%s/%s_definedRegion.pdf", outDir,g),
        dmrTrack = dmrGR,
        collapseTranscripts = "meta" # see Gviz documentation for GeneRegionTrack
    )
}



}, error = function(ex) {
    print(ex)
}, finally= {
    sink(NULL)
})

#' get range to plot methylation traces
#' @param df data.frame with DMRs and annotation
#' @param gene gene name to get the range for
#' @return (GRanges) range to plot
getPlottingRange <- function(df, gene){
    dmr <- df[which(df$Mannens_peakGene == gene),]
    if (nrow(dmr) == 0) {
        stop(sprintf("No DMRs found for gene %s", gene))  
    } else {
        gr <- GRanges(dmr$seqnames[1],
            IRanges(min(dmr$start), max(dmr$end))
        )
        return(gr)
    }
}
