# identify nearest gene to superenhancers
rm(list=ls())

source("../../utils_PaiLab.R")
require(GenomicRanges)

supEnhFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_FetalCBL_ChipSeq/CBL_Chipseq/Aldinger_FetalCBChipseq_superEnhancers.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_FetalCBL_ChipSeq/output/SE"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive=FALSE)
}

logFile <- sprintf("%s/log_annotateSuperEnhancers.txt", outDir)
sink(logFile, append=FALSE, split=TRUE)

tryCatch({
se <- read.delim(supEnhFile, sep="\t", header=FALSE)
cat(sprintf("Read %i super-enhancers\n", nrow(se)))

se_GR <- GRanges(seqnames = se$V1,
                 ranges = IRanges(start=se$V2, end=se$V3))

cat("Getting nearest genes\n")            
x <- getNearestGene(se_GR, gene_types="protein_coding",verbose=TRUE)

cat("Annotate these\n")
ndev <- get_neurodevGenes()
x$IsNeuroDev <- x$nearestGene %in% ndev
cat(sprintf("Found %i neurodev genes\n", sum(x$IsNeuroDev)))

g34 <- get_g34genes()
x$IsG34 <- x$nearestGene %in% g34
cat(sprintf("Found %i G34 genes\n", sum(x$IsG34)))


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