rm(list=ls())
library(BSgenome.Hsapiens.UCSC.hg38.masked) # needed for genNullSeqs
library(tidyr)
library(IRanges)
library(doParallel)
library(rtracklayer)
library(ggplot2)
library(Seurat)

source("../../utils_PaiLab.R")
source("utils.R")
source("getGRanges_OLenrichment.R")

cpg_dmr_date <- get_configs("CPG_DMR_DATE")
projectRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/"
projectRoot <- sprintf("%s/%s", projectRoot, cpg_dmr_date)
rlEPlinks <-  "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260407_Integrate/LinkPeaks_RLonly_output.RData"
geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

###rnaFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260402/Sarropoulos_RNA_Seurat.rds"

ciceroFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260406_CiceroRLonly/ciceroRLOnly_passCutoff.Rdata"
ciceroEPFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260410_DMR_overlap_EPs/Sarropoulos_RL_Cicero_inferredEP_upregGenes_oneTSS.txt"

numPerm <- 1000L
numCores <- 8L

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026"
outDir <- sprintf("%s/%s_DMROverlap", outDir, format(Sys.Date(),"%y%m%d"))
if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)

negDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRoverlap2/260401/negs"

logFile <- sprintf("%s/DMROverlap_Stats.log", outDir)
sink(logFile, append=FALSE, split=TRUE)

#' count the number of neurodev genes and MB genes in the overlap between tgtGR and ABC Enhancer-gene links
#' @param tgtGR GRanges of target regions (e.g. DMRs, negative controls)
#' @param rlGR GRanges of Sarropoulos enhancer-gene links, with gene name in a metadata column (e.g. "V4")
#' @param neurodevGenes vector of neurodevelopmental genes
#' @param mbGenes vector of medulloblastoma genes
#' @return (list) 1) numNeurodevOL: number of neurodev genes in the overlap, 
# 2) numMBOL: number of MB genes in the overlap, 3) numTotalOL: total number of genes in the overlap
countRL_geneStats <- function(tgtGR, rlGR, neurodevGenes, mbGenes) {
    combined <- GenomicRanges::findOverlaps(tgtGR, rlGR)
    tgtDF <- as.data.frame(tgtGR)
    rlDF <- as.data.frame(rlGR)
    x1 <- tgtDF[queryHits(combined),]
    x2 <- rlDF[subjectHits(combined),]
    both <- cbind(x1, x2)
    olGenes <- unique(both$gene_symbol) # assuming gene name is in column V4 of abcGR
    numNeurodevOL <- sum(olGenes %in% neurodevGenes)
    numMBOL <- sum(olGenes %in% mbGenes)
    numTotalOL <- length(olGenes)

    return(list(numNeurodevOL=numNeurodevOL, numMBOL=numMBOL, numTotalOL=numTotalOL))
}

tryCatch({

cat(sprintf("Reading DMRs from %s\n", cpg_dmr_date))
dmr <- getDMRs(projectRoot)

dmr <- regioneR::filterChromosomes(dmr, organism="hg", chr.type="canonical")
cat(sprintf("Filter standard chrom: DMRs left=%i\n", length(dmr)))

hg38 <- BSgenome.Hsapiens.UCSC.hg38.masked
hg38 <- keepStandardChroms(hg38) # remove alternate chroms 

cat("loading RL peak-gene links from Sarropoulos\n")
load(rlEPlinks)
rlGR <- LinkPeaks_output
cat(sprintf("Read %i RL peak-gene links\n", nrow(rlGR)))

# make a scatter plot of rlGR score and -log10(pval) to choose a threshold for filtering.
rlDF <- as.data.frame(rlGR)
p <- ggplot(as.data.frame(rlDF), aes(x=score, y=-log10(pvalue))) + geom_point() 
p <- p + xlab("RL score") + ylab("-log10(pvalue)") + ggtitle("Scatter plot of RL scores and p-values")
p <- p + geom_hline(yintercept=-log10(0.05), color="red", linetype="dashed") + 
    geom_vline(xintercept=0.3, color="blue", linetype="dashed")
p <- p + theme(axis.text=element_text(size=18))
ggsave(sprintf("%s/RL_score_pval_scatter.png", outDir), p, width=6, height=4)

genes <- rtracklayer::readGFF(geneFile)
genes <- subset(genes, gene_type %in% "protein_coding" & type == "gene")
rlGR <- rlGR[which(rlGR$gene_symbol %in% genes$gene_name)]
cat(sprintf("RL peak-gene links in protein-coding genes: %i\n", length(rlGR)))

# let's only look at the links for genes upregulated in the RL.
if (FALSE) { ### cuts statistical power too much.
    cat("Reading gene expression data\n")
    rna <- readRDS(rnaFile)
    DefaultAssay(rna) <- "SCT"
    Idents(rna) <- rna$precisest_label
    t0 <- Sys.time()
    rl_markers <- FindMarkers(
      object = rna,
      ident.1 = "progenitor_RL",
      ident.2 = NULL,
      min.pct = 0.1,
      test.use = 'wilcox'
    )
    print(Sys.time()-t0)
    cat(sprintf("Total genes tested for DEG = %s; Q < 0.1 = %s\n", 
        prettyNum(nrow(rl_markers), big.mark = ","), 
        prettyNum(sum(rl_markers$p_val_adj < 0.1), big.mark = ","))
    )
    upreg_genes <- rownames(rl_markers)[rl_markers$p_val_adj < 0.2 & rl_markers$avg_log2FC > 0]
    cat(sprintf("Total %s genes upregulated in RL progenitors\n", 
        prettyNum(length(upreg_genes), big.mark = ",")))

    cat("Filtering RL peak-gene links for those with genes upregulated in the RL\n")
    rlGR <- rlGR[which(rlGR$gene_symbol %in% upreg_genes)]
    cat(sprintf("RL peak-gene links with genes upregulated in the RL: %i\n", length(rlGR)))
}

###cat("----------------------------\n")
###cat("METHOD 1: LinkPeaks() from rhombic lip ATAC+RNA integration\n")
###cat("----------------------------\n")
###cat("DMR overlaps with RL E-P links\n")
###ol_EP <- getGRanges_OLenrichment(
###        pos=dmr,tgtGR=rlGR, numPerm=numPerm, negDir=negDir,
###        rngSeed=12345,genome=hg38, outDir=outDir,
###        tgtName="Sarropoulos RL peak-gene links", verbose=TRUE)
###
if (FALSE) { 
    cat("Reading genes of interest\n")
    neurodev_genes <- get_neurodevGenes()
    g34_genes <- get_g34genes()  

    cat("---------------------\n")
    cat("DMR links to neurodev genes\n")

    olreal <- countRL_geneStats(dmr, rlGR, neurodev_genes, g34_genes)
    negs <- dir(path=negDir, pattern="neg")
    negs <- negs[grep("(?<!unfiltered)\\.bed", negs, perl = TRUE)] # to use only filtered ones
    negs <- sprintf("%s/%s",negDir,negs)[1:250]
    cat(sprintf("Found %i neg sets\n",length(negs)))

    negStats <- data.frame(numNeurodevOL=integer(length(negs)), 
        numMBOL=integer(length(negs)), 
        numTotalOL=integer(length(negs)))

    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    nullol <- foreach (i=1:length(negs),.packages=c("GenomicRanges","rtracklayer")) %dopar% {  
        negGR <- import(negs[i])
        stats <- countRL_geneStats(negGR, rlGR, neurodev_genes, g34_genes)
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
}

cat("----------------------------\n")
cat("METHOD 2: Cicero co-accessibility peaks in rhombic lip\n")
cat("----------------------------\n")
load(ciceroFile)
conns <- cicero.good

###cat("\nCicero Peak 1\n")
###p1 <- do.call("rbind",strsplit(conns$Peak1, split="-"))
###p1 <- p1[!duplicated(p1),]
###gr1 <- makeGRangesFromDataFrame(data.frame(
###    "chr" = p1[,1],
###    "start" = as.numeric(p1[,2]),
###    "end" = as.numeric(p1[,3])
###))
###cat(sprintf("DMR overlaps with Cicero peak 1 (N=%i peaks)\n", length(gr1)))
###ol_p1 <- getGRanges_OLenrichment(
###        pos=dmr,tgtGR=gr1, numPerm=numPerm, negDir=negDir,
###        rngSeed=12345,genome=hg38, outDir=outDir,
###        tgtName="Sarropoulos RL Cicero Peak 1", verbose=TRUE)
###
###cat("\nCicero Peak 2\n")
###p2 <- do.call("rbind",strsplit(as.character(conns$Peak2), split="-"))
###p2 <- p2[!duplicated(p2),]
###gr2 <- makeGRangesFromDataFrame(data.frame(
###    "chr" = p2[,1],
###    "start" = as.numeric(p2[,2]),
###    "end" = as.numeric(p2[,3])
###))
###cat(sprintf("DMR overlaps with Cicero peak 2 (N=%i peaks)\n", length(gr2)))
###ol_p2 <- getGRanges_OLenrichment(
###        pos=dmr,tgtGR=gr2, numPerm=numPerm, negDir=negDir,
###        rngSeed=12345,genome=hg38, outDir=outDir,
###        tgtName="Sarropoulos RL Cicero Peak 2", verbose=TRUE)
###
cat("now reading the inferred EP pairs linked to")
ciceroEP <- read.delim(ciceroEPFile, header=TRUE)[,1:3]
ciceroEP <- ciceroEP[!duplicated(ciceroEP),]
gr3 <- GRanges(ciceroEP[,1], IRanges(ciceroEP[,2], ciceroEP[,3]))
ol_cep<- getGRanges_OLenrichment(
        pos=dmr,tgtGR=gr3, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="Sarropoulos RL Cicero EP", verbose=TRUE)

}, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
}, finally = {
    sink()
})