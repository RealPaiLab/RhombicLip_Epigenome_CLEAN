rm(list=ls())
library(BSgenome.Hsapiens.UCSC.hg38.masked) # needed for genNullSeqs
library(tidyr)
library(IRanges)
library(doParallel)

source("../../utils_PaiLab.R")
source("utils.R")

cpg_dmr_date <- get_configs("CPG_DMR_DATE")
projectRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/"
projectRoot <- sprintf("%s/%s", projectRoot, cpg_dmr_date)

dmrFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"
supEnhFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_FetalCBL_ChipSeq/CBL_Chipseq/Aldinger_FetalCBChipseq_superEnhancers.txt"
abcFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/Nasser-Neuronal-ABC_creTarget_hg38.bed"
matchedDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_CpGcompare/260602"

dt <- format(Sys.Date(),"%y%m%d")
if (!is.null(matchedDir)) {
    outDir <- sprintf("%s/DMRoverlap2_matched", projectRoot)
} else { 
    outDir <- sprintf("%s/DMRoverlap2",projectRoot)
}


if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)

negDir <- sprintf("%s/DMRoverlap2/260601/negs",projectRoot)
logFile <- sprintf("%s/DMROverlap_Stats.log", outDir)

numPerm <- 50L #1000L
numCores <- 10L
source("getGRanges_OLenrichment.R")

#' count the number of neurodev genes and MB genes in the overlap between tgtGR and ABC Enhancer-gene links
#' @param tgtGR GRanges of target regions (e.g. DMRs, negative controls)
#' @param abcGR GRanges of ABC enhancer-gene links, with gene name in a metadata column (e.g. "V4")
#' @param neurodevGenes vector of neurodevelopmental genes
#' @param mbGenes vector of medulloblastoma genes
#' @return (list) 1) numNeurodevOL: number of neurodev genes in the overlap, 
# 2) numMBOL: number of MB genes in the overlap, 3) numTotalOL: total number of genes in the overlap
countABC_geneStats <- function(tgtGR, abcGR, neurodevGenes, mbGenes) {
    combined <- GenomicRanges::findOverlaps(tgtGR, abcGR)
    tgtDF <- as.data.frame(tgtGR)
    abcDF <- as.data.frame(abcGR)
    x1 <- tgtDF[queryHits(combined),]
    x2 <- abcDF[subjectHits(combined),]
    both <- cbind(x1, x2)
    olGenes <- unique(both$V4) # assuming gene name is in column V4 of abcGR
    numNeurodevOL <- sum(olGenes %in% neurodevGenes)
    numMBOL <- sum(olGenes %in% mbGenes)
    numTotalOL <- length(olGenes)

    return(list(numNeurodevOL=numNeurodevOL, numMBOL=numMBOL, numTotalOL=numTotalOL))
}

logFile <- sprintf("%s/DMRoverlap_Stats2.log",outDir)
sink(logFile, append=FALSE, split=TRUE)

tryCatch({

cat(" *** Matched negative indices provided\n *** ")
cat(sprintf("Reading DMRs from %s\n", cpg_dmr_date))
dmrs <- read.delim(dmrFile, header=T, stringsAsFactors = FALSE)
dmrs$dmr_name <- sprintf("%s-%s-%s", dmrs$chr, dmrs$start, dmrs$end)
cat(sprintf("Read %i DMRs\n", nrow(dmrs)))

dmr_df <- dmrs
dmr <- GRanges(seqnames = dmrs$chr,
                 ranges = IRanges(start=dmrs$start, end=dmrs$end),
                 dmr_name = dmrs$dmr_name)

#dmr <- regioneR::filterChromosomes(dmr, organism="hg", chr.type="canonical")
cat(sprintf("Filter standard chrom: DMRs left=%i\n", length(dmr)))

hg38 <- BSgenome.Hsapiens.UCSC.hg38.masked
hg38 <- keepStandardChroms(hg38) # remove alternate chroms 

cat("HARs\n")
cat("----------------\n")
hars <- getHARs()
ol <- findOverlaps(dmr, hars)
dmr_ol <- dmr_df[which(dmr_df$dmr_name %in% dmr$dmr_name[queryHits(ol)]),]
cat(sprintf("Num DMRs overlapping HARs: %i\n", length(unique(queryHits(ol)))))
cat(sprintf("# hypoDMRs = %i, # hyperDMRs = %i\n", sum(dmr_ol$diff.Methy<0), sum(dmr_ol$diff.Methy>0)))
ol_sv <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=hars, numPerm=numPerm, negDir=negDir, idxDir=matchedDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="HARs", verbose=TRUE)

cat("Super-enhancers\n")
cat("----------------\n")
se <- read.delim(supEnhFile, sep="\t", header=FALSE)
cat(sprintf("Read %i super-enhancers\n", nrow(se)))

se_GR <- GRanges(seqnames = se$V1,
                 ranges = IRanges(start=se$V2, end=se$V3))       

ol <- findOverlaps(dmr, se_GR)
dmr_ol <- dmr_df[which(dmr_df$dmr_name %in% dmr$dmr_name[queryHits(ol)]),]
cat(sprintf("Num DMRs overlapping super-enhancers: %i\n", length(unique(queryHits(ol)))))
cat(sprintf("# hypoDMRs = %i, # hyperDMRs = %i\n", sum(dmr_ol$diff.Methy<0), sum(dmr_ol$diff.Methy>0)))
superenh_OL <- cbind(as.data.frame(dmr[queryHits(ol)]), as.data.frame(se_GR[subjectHits(ol)]))
superenh_OL <- superenh_OL[!duplicated(superenh_OL),]
write.table(superenh_OL, file=sprintf("%s/DMR_superEnhancer_overlap.txt", outDir), sep="\t", quote=F, row.names=F)

ol_sv <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=se_GR, numPerm=numPerm, negDir=negDir, idxDir=matchedDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="SuperEnhancers", verbose=TRUE)

# for ABC genes we need to compute overlap for each negative set separately, and count 
# the number of neurodev genes & MB genes.
cat("Reading ABC genes\n")
abc <- read.delim(abcFile, header=FALSE, stringsAsFactors = FALSE)
abc <-  subset(abc, V1 %in% paste("chr", c(1:22, "X", "Y"),sep=""))
abcGR <- makeGRangesFromDataFrame(
    abc, 
    keep.extra.columns = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3"
)

cat("Reading genes of interest\n")
neurodev_genes <- get_neurodevGenes()
g34_genes <- get_g34genes()  

olreal <- countABC_geneStats(dmr, abcGR, neurodev_genes, g34_genes)
negs <- dir(path=negDir, pattern="neg")
negs <- negs[grep("(?<!unfiltered)\\.bed", negs, perl = TRUE)] # to use only filtered ones
negs <- sprintf("%s/%s",negDir,negs)
cat(sprintf("Found %i neg sets\n",length(negs)))

negStats <- data.frame(numNeurodevOL=integer(length(negs)), 
    numMBOL=integer(length(negs)), 
    numTotalOL=integer(length(negs)))

cl <- makeCluster(numCores)
registerDoParallel(cl)
on.exit(stopCluster(cl))

nullol <- foreach (i=1:length(negs),.packages=c("GenomicRanges","rtracklayer")) %dopar% {  
    negGR <- import(negs[i])
    stats <- countABC_geneStats(negGR, abcGR, neurodev_genes, g34_genes)
    negStats[i,] <- unlist(stats)
    #cat(sprintf("Neg %i/%i: %s - numNeurodevOL=%i, numMBOL=%i, numTotalOL=%i\n", i, length(negs), basename(negs[i]), stats$numNeurodevOL, stats$numMBOL, stats$numTotalOL))
}
negStats <- do.call("rbind", nullol)

# first test for neurodev genes.
# strict one tail greater-than test
cat("Neurodev gene overlap\n")
nullol <- negStats[,"numNeurodevOL"]
realol <- olreal$numNeurodevOL
numPerm <- length(nullol)
tgtName <- "NeurodevGenes"
pval <- max(sum(nullol>=realol)/numPerm, 1/numPerm)
cat(sprintf("Num. real = %i, Mean null = %1.2f\n", realol, mean(nullol)))
cat(sprintf("Num null >= real overlaps = %i (p < %1.2e)\n", 
    sum(nullol>=realol), pval))
perms <- data.frame(iter=1:numPerm, nullol=nullol)
p <- ggplot(perms, aes(x=nullol))+ geom_density() +
    geom_vline(xintercept=realol,color="red", size = 1) +
    xlab("Percent ranges overlap") + 
    ggtitle(sprintf("Num. mapping to %s (N=%i permutations)", tgtName, numPerm)) +
    theme(axis.text=element_text(size=18))
ggsave(sprintf("%s/%s_neurodev_overlap.png", outDir, tgtName), p, width=6, height=4)
  

## test for MB genes
cat("MB gene overlap\n")
nullol <- negStats[,"numMBOL"]
realol <- olreal$numMBOL
numPerm <- length(nullol)
pval <- max(sum(nullol>=realol)/numPerm, 1/numPerm)
tgtName <- "MBGenes"
cat(sprintf("Num. real = %i, Mean null = %1.2f\n", realol, mean(nullol)))
cat(sprintf("Num null >= real overlaps = %i (p < %1.2e)\n", 
    sum(nullol>=realol), pval)) 
p <- ggplot(perms, aes(x=nullol))+ geom_density() +
    geom_vline(xintercept=realol,color="red", size = 1) +
    xlab("Percent ranges overlap") + 
    ggtitle(sprintf("Num. mapping to %s (N=%i permutations)", tgtName, numPerm)) +
    theme(axis.text=element_text(size=18))
ggsave(sprintf("%s/%s_mb_overlap.png", outDir, tgtName), p, width=6, height=4)

}, error=function(ex){
    print(ex)
}, finally={
    cat("Done.\n")
    sink(NULL)
})

