#  Identify TF-enhancer-target gene triplets using predicted TFBS overlapping DMRs

rm(list=ls())

source("../../utils_PaiLab.R")
require(GenomicRanges)

AMEdir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711"
predictedTFBS1 <- sprintf("%s/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/sequences.tsv",AMEdir)
predictedTFBS2 <- sprintf("%s/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/sequences.tsv",AMEdir)

ameF1 <- sprintf("%s/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv",AMEdir)
ameF2 <- sprintf("%s/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv",AMEdir)

inDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_link2Genes_ABC"

inFile <- sprintf("%s/250729/DMR_link2Genes_ABC_250729.sif1_1.sif",inDir) 

fullDMR2geneFile <- sprintf("%s/250729/DMR_AnnotatedAll_250729.tsv",inDir)

dt <- format(Sys.Date(),"%y%m%d")
outDir <- dirname(inFile)
outFile <- sprintf("%s/DMR_link2Genes_ABC_annotatedWithTFs_%s.csv",outDir,dt)

logFile <- sprintf("%s/DMR_link2Genes_AnnotateDMRs_%s.log",outDir,dt)
cat(sprintf("Logging to %s\n",logFile))

sink(logFile, split=TRUE)

tryCatch({

dmr2gene <- read.delim(inFile,sep="\t",header=F)
colnames(dmr2gene) <- c("target","type","source")
cat(sprintf("Read %i DMR-gene links\n",nrow(dmr2gene)))

crsplit <- do.call(rbind, strsplit(dmr2gene$source,"[:-]"))
enhGR <- makeGRangesFromDataFrame(
    data.frame(
        chr = crsplit[,1],  
        start = as.integer(crsplit[,2]),
        end = as.integer(crsplit[,3]),
        name = dmr2gene$source
    ), keep.extra.columns = TRUE)


cat("\n\nProcessing TF motifs from AME results\n")
tfbs1 <- read.delim(predictedTFBS1, header=TRUE, stringsAsFactors = FALSE)
tfbs2 <- read.delim(predictedTFBS2, header=TRUE, stringsAsFactors = FALSE)
tfbs <- rbind(tfbs1, tfbs2)
tfbs <- subset(tfbs, class %in% "tp")[,c(2,4)]
cat(sprintf("Read %i predicted TFBS\n", nrow(tfbs)))
tfbs <- tfbs[!duplicated(tfbs),]
cat(sprintf("%i unique TF motifs\n", length(unique(tfbs$motif_ID))))
cat(sprintf("After removing duplicates, %i predicted TFBS\n", nrow(tfbs)))
tf <- unique(tfbs$motif_ID)
tfbs[,1] <- substr(tfbs[,1], 1, regexpr("\\.", tfbs[,1])-1)
cat(sprintf("After removing sub-motif info, %i unique TF motifs\n", length(unique(tfbs$motif_ID))))
tfbs <- tfbs[!duplicated(tfbs),]
cat(sprintf("After removing duplicate motif IDs, %i predicted TFBS\n", nrow(tfbs)))
tfbs[,2] <- sub(":", "-", tfbs[,2])
#split column 2 by "-"
tfbs_split <- do.call(rbind, strsplit(tfbs[,2], "-"))
tfbs_gr <- makeGRangesFromDataFrame(
    data.frame(
        chr = tfbs_split[,1],  
        start = as.integer(tfbs_split[,2]),
        end = as.integer(tfbs_split[,3]),
        motif_ID = tfbs[,1]
    ), keep.extra.columns = TRUE)

# find overlaps between enhGR and tfbs_gr
ol <- findOverlaps(enhGR, tfbs_gr)
enhGR_tfbs <- enhGR[queryHits(ol)]
enhGR_tfbs$motif_ID <- tfbs_gr[subjectHits(ol)]$motif_ID
enhGR_tfbs <- as.data.frame(enhGR_tfbs)
cat(sprintf("Found %i overlaps between enhancers and predicted TFBS\n", nrow(enhGR_tfbs)))
enhGR_tfbs <- enhGR_tfbs[!duplicated(enhGR_tfbs[,c("name","motif_ID")]),]
cat(sprintf("After removing duplicates, %i unique enhancer-TFBS pairs\n", nrow(enhGR_tfbs)))


cat("\n\nFiltering by genes DEG in the rhombic lip\n")
RL_DEG <- get_configs("RL_DEG")
xprFull <- read.table(RL_DEG, header=TRUE, stringsAsFactors = FALSE)
xpr <- subset(xprFull, FDR < 0.05)
cat(sprintf("Read %i DEGs with PValue < 0.05\n", nrow(xpr)))
xpr$geneName <- rownames(xpr)

enhGR_tfbs_deg <- merge(enhGR_tfbs, xpr, by.x="motif_ID", by.y="geneName")
cat(sprintf("\n\nAfter merging with DEGs, %i enhancer-TFBS-DEG rows\n", nrow(enhGR_tfbs_deg)))
# how many unique motif_IDs?
cat(sprintf("There are %i unique motif_IDs in the merged data\n", 
length(unique(enhGR_tfbs_deg$motif_ID))))
# how many unique DMR locations (name column)?
cat(sprintf("There are %i unique DMR locations in the merged data\n", 
length(unique(enhGR_tfbs_deg$name))))

enhGR_tfbs_deg <- enhGR_tfbs_deg[c("name","motif_ID")]
colnames(enhGR_tfbs_deg) <- c("source","motif_ID")

x <- merge(dmr2gene, enhGR_tfbs_deg, by.x="source", by.y="source",
    all.x=TRUE)
write.table(x, file=outFile, sep=",", row.names=FALSE, quote=FALSE )

x2 <- x[,c("target","type","motif_ID")]
colnames(x2) <- c("TargetGene","InteractionType","TF")
x2 <- x2[!duplicated(x2),]
outFile2 <- sub(".csv","_uniqueInteractions.txt", outFile)
write.table(x2, file=outFile2, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE )
cat(sprintf("Wrote %i unique gene-DMR-TF interactions to %s\n", nrow(x2), outFile2))
cat("Done.\n")

###source("../../../../../MISC/JASPAR_getTFBS/getTFBSMotifMatrix.R")
###  jaspar_hg38 <- get_configs("JASPAR_DB")
###
###curwd <- getwd()
###  
###setwd("../../../../../MISC/JASPAR_getTFBS/")
###  tfbs_file <- sprintf("%s/dmr_jaspar_%s.out", outDir, dt)
###  if (file.exists(tfbs_file)) file.remove(tfbs_file)
###  jaspar_res <- getTFBSMotifMatrix(enhGR, 
###                     jaspar=jaspar_hg38, 
###                     outFile = tfbs_file, 
###                     tmpDir = outDir, 
###                     convertMat = TRUE
###  )
###  setwd(curwd) 


cat("\n\nNow annotate TFs for the full set of DMR-gene links\n")
fullDMR2gene <- read.delim(fullDMR2geneFile, 
    sep="\t", header=TRUE, stringsAsFactors = FALSE
)
cat(sprintf("Read %i full DMR-gene links\n", nrow(fullDMR2gene)))

dmrcoord <- do.call(rbind, strsplit(fullDMR2gene$name,"[:-]"))
fullEnhGR <- makeGRangesFromDataFrame(
    data.frame(
        chr = dmrcoord[,1],  
        start = as.integer(dmrcoord[,2]),
        end = as.integer(dmrcoord[,3]),
        name = fullDMR2gene$name
    ), keep.extra.columns = TRUE)
ol2 <- findOverlaps(fullEnhGR, tfbs_gr)
fullEnhGR_tfbs <- fullEnhGR[queryHits(ol2)]
fullEnhGR_tfbs$motif_ID <- tfbs_gr[subjectHits(ol2)]$motif_ID
fullEnhGR_tfbs <- as.data.frame(fullEnhGR_tfbs)
cat(sprintf("Found %i overlaps between full enhancers and predicted TFBS\n", nrow(fullEnhGR_tfbs)))
fullEnhGR_tfbs <- fullEnhGR_tfbs[!duplicated(fullEnhGR_tfbs[,c("name","motif_ID")]),]
cat(sprintf("After removing duplicates, %i unique enhancer-TFBS pairs\n", nrow(fullEnhGR_tfbs)))

xpr2 <- subset(xprFull, PValue < 0.05)
xpr2$geneName <- rownames(xpr2)
cat("For the full set, using DEG cutoff of PValue < 0.05\n")
fullEnhGR_tfbs_deg <- merge(fullEnhGR_tfbs, xpr2, by.x="motif_ID", by.y="geneName")
cat(sprintf("\n\nAfter merging with DEGs, %i full enhancer-TFBS-DEG rows\n", nrow(fullEnhGR_tfbs_deg)))
# how many unique motif_IDs?
cat(sprintf("There are %i unique motif_IDs in the merged full data\n", 
length(unique(fullEnhGR_tfbs_deg$motif_ID))))
# how many unique DMR locations (name column)?
cat(sprintf("There are %i unique DMR locations in the merged full data\n", 
length(unique(fullEnhGR_tfbs_deg$name))))   
fullEnhGR_tfbs_deg <- fullEnhGR_tfbs_deg[c("name","motif_ID")]
colnames(fullEnhGR_tfbs_deg) <- c("source","motif_ID")
x_full <- merge(fullDMR2gene, fullEnhGR_tfbs_deg, by.x="name", by.y="source",
    all.x=TRUE)
outFile_full <- sub(".csv","_fullSetOfDMRs2genes.csv", outFile)
cat(sprintf("\n\nFull set of TF-DMR-gene links has %i unique DMRs, %i unique TF motifs, and %i unique target genes.\n Total number of links = %i\n",
    length(unique(x_full$name)),
    length(unique(x_full$motif_ID)),
    length(unique(x_full$ABC_gene)),
    nrow(x_full)))
write.table(x_full, file=outFile_full, sep=",", row.names=FALSE, quote=FALSE)


},error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")
},finally={
    sink()
})