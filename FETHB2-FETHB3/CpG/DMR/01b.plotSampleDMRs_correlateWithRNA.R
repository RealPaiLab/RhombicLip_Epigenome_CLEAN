# plot detailed view of selected DMRs and correlate promoter methylation with gene expression 
rm(list=ls())
library(ggplot2)
library(bsseq)
library(DSS)
require(reshape2)
library(EnvStats)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("locusPlotter.R")
source("../../utils_PaiLab.R")
source("getM_AllSamples_GRanges.R")

cytoDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/CpG_snpFiltered"
CPG_DMR_FILE <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"

phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv"

rnaDEGfile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/231123/edgeR_RLVZvsSVZ_231123.txt"

olFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs_cCRE_overlap_240712.txt"


outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_plots"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/plotSampleDMRs.log", outDir)
sink(logFile, split=TRUE)
tryCatch({
  

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

cat("* Reading DMRs\n")
dmrs <- read.delim(
    CPG_DMR_FILE, 
    header = TRUE, sep="\t",
    stringsAsFactors = FALSE
)

dmrs$name <- paste(dmrs$chr, 
                      dmrs$start, 
                      dmrs$end, 
                      sep = "-")
#dmrs <- dmrs[order(dmrs$diff.Methy, decreasing = TRUE),]

ol <- read.delim(olFile, 
    header = TRUE, sep="\t",
    stringsAsFactors = FALSE
)

tss_gr <- get_gencode_anno(gene_type="protein_coding")
tss_gr <- promoters(tss_gr, 
    upstream = 1000, 
    downstream = 150
)

cat("Getting promoter-level DNA methylation for all samples\n")
t0 <- Sys.time()
prom_meth <- getM_AllSamples_GRanges(
  cytoDir, 
  tss_gr
)
colnames(prom_meth$pctM) <- tss_gr$name
print(Sys.time() - t0)
save(prom_meth,tss_gr,
     file = sprintf("%s/promoter_methylation.RData", outDir))
pctM <- as.data.frame(prom_meth$pctM)


cat("* Analyzing methylation in promoter regions overlapping DMRs\n")
ol_gr <- makeGRangesFromDataFrame(
  ol, 
  keep.extra.columns = TRUE, 
  seqnames.field = "DMR.seqnames", 
  start.field = "DMR.start", 
  end.field = "DMR.end"
)
names(ol_gr) <- paste(ol$DMR.seqnames, 
                      ol$DMR.start, 
                      ol$DMR.end, 
                      sep = "-")                                    
ol$name <- names(ol_gr)                      

# create a table of DMRs with their overlap with TSS
dmr_tss <- findOverlaps(ol_gr, tss_gr)
both <- data.frame(
  dmr = names(ol_gr)[queryHits(dmr_tss)],
  tss = tss_gr$name[subjectHits(dmr_tss)],
  stringsAsFactors = FALSE
  )
ol2 <- merge(ol, both, by.x = "name", by.y = "dmr", all.x = TRUE)
ol2 <- merge(x=ol2, y=dmrs, 
              by="name", all.x = TRUE)
ol2 <- ol2[!duplicated(ol2),]                 
prom <- subset(ol2, ol2$ENCODE.cCRE.Type %in% c("PLS","PLS,CTCF-bound") 
        & !is.na(ol2$tss))
prom <- prom[!duplicated(prom),]      

# ATOH7, BRINP1, CRMP1, NHLH1, NTRK2, STMN1, USP9X, SPTBN2, FOXJ1, GPR12

# define colours for plots
if (all.equal(pheno$sample, sampleNames(bsObj$bs))!=TRUE) {
  stop("Sample names in pheno and bsObj do not match!")
}
cols <- c("#7B3294","#008837")[factor(pheno$ROI,
  levels=c("VZ","SVZ"))]

x2 <- prom[,c("name","tss","diff.Methy")]
x2 <- x2[!duplicated(x2),]
comb <- 
  merge(
    x = x2,
    y = DEG[,c("gene","logFC","PValue","FDR")],
    by.x = "tss", 
    by.y = "gene"
  )

DEG <- read.delim(
  rnaDEGfile, 
  header = TRUE, 
cat("* Reading DEG results ...\n")
  stringsAsFactors = FALSE

)
DEG$gene <- rownames(DEG)
DEG <- DEG[order(DEG$FDR),]

# plot a scatterplot of diffMethy vs logFC, add a fitted line and text indicating correlation and significance
p <- ggplot(comb, aes(x = logFC, y = diff.Methy*100)) +
  geom_point() + geom_smooth(method="lm") +
  labs(title = "Diff Methylation vs Log2 Fold Change",
       y = "% DNA methylation increase in RL-VZ", 
       x = "Log2 Fold Change in RL-VZ") +
  theme_minimal(base_size = 12)
p <- p + stat_cor(method="spearman",label.y=25, label.x=1)
p <- p + geom_text_repel(
  data = subset(comb, abs(diff.Methy) > 0.05 & abs(logFC) > 0.3),
  aes(label = tss),
  size = 3,
  nudge_x = 0.1,
  nudge_y = 0.1,
  show.legend = FALSE
)
ggsave(
  filename = sprintf("%s/DMR_vs_DEG_promoters.pdf", outDir),
  plot = p,
  width = 8, height = 6
)

# clean pctM
pctM <- as.data.frame(prom_meth$pctM)
dpos <- regexpr("-",rownames(pctM))
rownames(pctM) <- substr(rownames(pctM), 1, dpos-1)
pctM$sample <- rownames(pctM)
pctM <- melt(pctM)

pctM <- merge(pctM, pheno[,c("sample","ROI")], 
  by="sample", all.x=TRUE)

# compute median difference VZ - SVZ for each gene
pctM$ROI <- factor(pctM$ROI, levels=c("VZ","SVZ"))
pctM_diff <- pctM %>% 
  group_by(variable, ROI) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  summarise(
    diffMethy = median_value[ROI == "VZ"] - median_value[ROI == "SVZ"]
  )

comb <- merge(x=DEG, y=pctM_diff, 
              by.x = "gene", by.y = "variable")
sig <- subset(comb, FDR < 0.05)
# plot a scatterplot of diffMethy vs logFC, add a fitted line and text indicating correlation and significance
p <- ggplot(sig, aes(x = logFC, y = diffMethy*100)) +
  geom_point() + geom_smooth(method="lm") +
  labs(title = "Diff Methylation vs Log2 Fold Change",
       y = "% DNA methylation increase in RL-VZ",
       x = "Log2 Fold Change in RL-VZ") +
  theme_minimal(base_size = 12) 
p <- p + stat_cor(method="pearson",label.y=25, label.x=1)
p <- p + geom_text_repel(
  data = subset(sig, abs(diffMethy) > 0.05 & abs(logFC) > 0.3),
  aes(label = gene),
  size = 3,
  nudge_x = 0.1,
  nudge_y = 0.1,
  show.legend = FALSE
)
ggsave(
  filename = sprintf("%s/AllPromoters_diffMethy_vs_Log2FC_sig.pdf",
    outDir),
  plot = p,
  width = 8, height = 6
)

cat("* Plot methylation for promoter regions of selected genes\n")
inFile <- get_configs("GENCODE_GENE_FILE") 
genes <- rtracklayer::readGFF(inFile)

# plot detailed methylation view for selected genes
geneList <- c("SOX2","WLS","EOMES", "NEUROD1","NEUROD2", "NHLH1",
  "CRB","TMSB10")
ups <- c(2000, 5000, 10000, 10000, 10000, 5000, 5000, 10000)
downs <- c(10000, 20000,20000, 40000, 40000, 10000, 5000, 10000)
for (i in c(1,6)) { # ,2)) { #1:length(geneList)) {
  g <- geneList[i]
  print(g)
  rg <- tss_gr[which(tss_gr$name == g)]
  rg <- promoters(rg, 
    upstream = ups[i], 
    downstream = downs[i])
  curg <- subset(genes, gene_name == g & type == "gene")
  curg$exon_number <- 1
  curg$isoforms <- 1

  methyPlotter(
    bsObj$bs, 
    region = rg, 
    extend = 0,
    geneTrack = curg,
    cols = cols, 
    lgd_cols = c("VZ" = "#7B3294", "SVZ" = "#008837"),
    geneName = g,
    outFile = sprintf("%s/%s.pdf", outDir,g)
  )
}

}, error = function(e) {
  cat("Error in setting options: ", e$message, "\n")
})







