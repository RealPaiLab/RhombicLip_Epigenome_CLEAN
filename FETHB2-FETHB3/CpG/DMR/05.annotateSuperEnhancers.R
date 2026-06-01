# identify nearest gene to superenhancers
rm(list=ls())

source("../../utils_PaiLab.R")
source("../overlapEnrichment/getGRanges_OLenrichment.R")
require(GenomicRanges)
require(rtracklayer)

supEnhFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_FetalCBL_ChipSeq/CBL_Chipseq/Aldinger_FetalCBChipseq_superEnhancers.txt"
outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_FetalCBL_ChipSeq/output/SE"
geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"


#' count the number of neurodev genes and MB genes in the overlap between tgtGR and ABC Enhancer-gene links
#' @param tgtGR GRanges of target regions (e.g. DMRs, negative controls)
#' @param geneGR GRanges of gene annotations, with gene name in a metadata column (e.g. "gene_name")
#' @param neurodevGenes vector of neurodevelopmental genes
#' @param mbGenes vector of medulloblastoma genes
#' @return (list) 1) numNeurodevOL: number of neurodev genes in the overlap, 
# 2) numMBOL: number of MB genes in the overlap, 3) numTotalOL: total number of genes in the overlap
count_geneStats <- function(tgtGR, geneGR, neurodevGenes, mbGenes) {
    
    n <- nearest(tgtGR, geneGR)
    tgtGR$nearestGene <- geneGR$gene_name[n] 

    tgtGR$IsNeuroDev <- tgtGR$nearestGene %in% ndev
    tgtGR$IsG34 <- tgtGR$nearestGene %in% g34
    return(list(numNeurodevOL=sum(tgtGR$IsNeuroDev), numMBOL=sum(tgtGR$IsG34), numTotalOL=length(unique(tgtGR$nearestGene))))
}

negDir <- sprintf("%s/260415/negs", outDir)

dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive=FALSE)
}


logFile <- sprintf("%s/log_annotateSuperEnhancers.txt", outDir)
numCores <- 10L

sink(logFile, append=FALSE, split=TRUE)

tryCatch({
se <- read.delim(supEnhFile, sep="\t", header=FALSE)
cat(sprintf("Read %i super-enhancers\n", nrow(se)))

se_GR <- GRanges(seqnames = se$V1,
                 ranges = IRanges(start=se$V2, end=se$V3))

cat("Getting nearest genes\n")            
x <- getNearestGene(se_GR, gene_types="protein_coding",verbose=TRUE)

numPerm <- 1000

ndev <- get_neurodevGenes()
g34 <- get_g34genes()

cat("Generating null sets\n")
t0 <- Sys.time()
if (!dir.exists(negDir)) {
  dir.create(negDir, recursive=FALSE)
  x <- getNullGRanges(se_GR, numPerm, outDir=negDir, numCores=10L, filterBismapMappable=TRUE)
}
print(Sys.time() - t0)

cat("Reading gene definitions\n")
genes <- rtracklayer::import(geneFile)
genes <- subset(genes, gene_type == "protein_coding" & type == "gene") 
genes$TSS <- genes$start
genes$TSS[which(genes$strand=="-")] <- genes$end[which(genes$strand=="-")]

olreal <- count_geneStats(se_GR, genes, neurodevGenes=ndev, mbGenes=g34)
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
    stats <- count_geneStats(negGR, genes, ndev, g34)
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
cat(sprintf("Num. real = %i, Median null = %1.2f\n", realol, median(nullol)))
cat(sprintf("Num null >= real overlaps = %i (p < %1.2e)\n", 
    sum(nullol>=realol), pval))


## test for MB genes
cat("MB gene overlap\n")
nullol <- negStats[,"numMBOL"]
realol <- olreal$numMBOL
numPerm <- length(nullol)
pval <- max(sum(nullol>=realol)/numPerm, 1/numPerm)
tgtName <- "MBGenes"
cat(sprintf("Num. real = %i, median null = %1.2f\n", realol, median(nullol)))
cat(sprintf("Num null >= real overlaps = %i (p < %1.2e)\n", 
    sum(nullol>=realol), pval)) 

# create dataframe with nearest gene and write to file
outFile <- sprintf("%s/annotated_superEnhancers.txt", outDir)
write.table(x, file=outFile, sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Wrote annotated super-enhancers to %s\n", outFile))

},error=function(ex){
    print(ex)
},
finally={
    sink()
})