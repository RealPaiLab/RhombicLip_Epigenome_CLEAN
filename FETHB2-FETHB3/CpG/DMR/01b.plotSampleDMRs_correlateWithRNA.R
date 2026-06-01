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

#CPG_DMR_FILE <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/moreStringent_251210/DMRs_moreStringent.csv"
#olFile <- sprintf("%s/DMRs_cCRE_overlap_251210.txt",dirname(CPG_DMR_FILE))

phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv"

rnaDEGfile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/251210/edgeR_RLVZvsSVZ_OlderThan14PCW_251210.txt"

olFile <- sprintf("%s/DMRs_cCRE_overlap_240712.txt",dirname(CPG_DMR_FILE))

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_plots"

dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

upstream <- 2000
downstream <- 1000
minCvg <- 5L
minSamp <- 9L

promFile <- sprintf("%s/promoter_methylation_up%i_down%i.RData", 
  outDir, upstream, downstream)

logFile <- sprintf("%s/plotSampleDMRs_minCvg%i_up%i_down%i.log", outDir, minCvg, upstream, downstream)
sink(logFile, split=TRUE)
tryCatch({

cat("Parameters:\n")
cat(sprintf("  DMR file: %s\n", CPG_DMR_FILE))
cat(sprintf("  Overlap file: %s\n", olFile))
cat(sprintf("  RNA DEG file: %s\n", rnaDEGfile))
cat(sprintf("  minCvg: %i\n", minCvg))
cat(sprintf("  minSamp: %i\n", minSamp))
cat(sprintf("  Promoter region: upstream %i, downstream %i\n", upstream, downstream))
cat("------------------------------------------------------------\n")

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

ol <- read.delim(olFile, 
    header = TRUE, sep="\t",
    stringsAsFactors = FALSE
)

if (!file.exists(promFile)) {

cat("* Reading BSseq object\n")
t0 <- Sys.time()
bsObj <- readBS(inDir = cytoDir, 
                  grepPattern = "\\.txt\\.gz$")
print(Sys.time() - t0)

tss_gr <- get_gencode_anno(gene_type="protein_coding")
tss_gr <- promoters(tss_gr, 
    upstream = upstream, 
    downstream = downstream
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
     file = promFile)
} else {
  cat("* Loading promoter-level DNA methylation for all samples\n")
  load(promFile)
}

cvg <- prom_meth$COV
idx <- which(colSums(cvg >= minCvg, na.rm=TRUE) >= minSamp)
cat(sprintf("* Retaining %i of %i promoters with coverage >= %i in at least %i samples\n", 
  length(idx), ncol(cvg), minCvg, minSamp))
prom_meth$pctM <- prom_meth$pctM[,idx]
prom_meth$COV <- prom_meth$COV[,idx]

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
###if (all.equal(pheno$sample, sampleNames(bsObj$bs))!=TRUE) {
###  stop("Sample names in pheno and bsObj do not match!")
###}
###cols <- c("#7B3294","#008837")[factor(pheno$ROI,
###  levels=c("VZ","SVZ"))]

x2 <- prom[,c("name","tss","diff.Methy")]
x2 <- x2[!duplicated(x2),]

cat("* Reading DEG results ...\n")
DEG <- read.delim(
  rnaDEGfile, 
  header = TRUE, 
  stringsAsFactors = FALSE

)
DEG$gene <- rownames(DEG)
DEG <- DEG[order(DEG$FDR),]

comb <- 
  merge(
    x = x2,
    y = DEG[,c("gene","logFC","PValue","FDR")],
    by.x = "tss", 
    by.y = "gene"
  )

# plot a scatterplot of diffMethy vs logFC, add a fitted line and text indicating correlation and significance
p <- ggplot(comb, aes(x = logFC, y = diff.Methy*100)) +
  geom_point() + geom_smooth(method="lm") +
  labs(title = sprintf("Diff Methylation vs Log2 Fold Change (Up: %1.1fkb, Down: %1.1fkb)\n(N=%i, MinCVG=%i, MinSamp=%i)", upstream/1000, downstream/1000, nrow(comb), minCvg, minSamp),
       y = "% DNA methylation increase in RL-VZ", 
       x = "Log2 Fold Change in RL-VZ") +
  theme_minimal(base_size = 24)
p <- p + stat_cor(method="pearson",label.y=25, label.x=1, size=6)
p <- p + geom_text_repel(
  data = subset(comb, abs(diff.Methy) > 0.05 & abs(logFC) > 0.3),
  aes(label = tss),
  size = 6,
  nudge_x = 0.1,
  nudge_y = 0.1,
  show.legend = FALSE
) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_vline(xintercept = 0, linetype="dashed", color = "red") 

colnames(comb)[2:3] <- paste0("DMR_", colnames(comb)[2:3])
colnames(comb)[4:6] <- paste("RL_DEG_", colnames(comb)[4:6], sep="")
write.table(comb,file=sprintf("%s/DMR_DEG_minCvg%i_up%i_down%i.txt", outDir, minCvg, upstream, downstream), 
  sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
browser()

cr <- cor.test(comb$logFC, comb$diff.Methy, method="pearson")
cat(sprintf("DMR vs DEG: Correlation between promoter methylation change and gene expression change:\n"))
cat(sprintf("  Number of promoters analyzed: %i\n", nrow(comb)))
cat(sprintf("  Pearson's r = %1.3f, p-value = %1.3e\n", 
  cr$estimate, cr$p.value)) 

tmp <- data.frame(
  minCvg=minCvg,
  minSamp=minSamp,
  upstream=upstream,
  downstream=downstream,
  n_promoters=nrow(comb),
  cor_pearson=cr$estimate,
  pval_pearson=cr$p.value,
  cor_spearman=cor.test(comb$logFC, comb$diff.Methy, method="spearman")$estimate,
  pval_spearman=cor.test(comb$logFC, comb$diff.Methy, method="spearman")$p.value
)

write.table(
  tmp,
  file = sprintf("%s/DMR_vs_DEG_promoters_correlation_summary_minCvg%i_up%i_down%i.txt", 
    outDir, minCvg, upstream,downstream),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

ggsave(
  filename = sprintf("%s/DMR_vs_DEG_promoters_minCvg%i_up%i_down%i.pdf", outDir, minCvg, upstream, downstream),
  plot = p,
  width = 8, height = 6
)

# clean pctM
pctM <- as.data.frame(prom_meth$pctM)
dpos <- regexpr("-",rownames(pctM))
rownames(pctM) <- substr(rownames(pctM), 1, dpos-1)
pctM$sample <- rownames(pctM)
pctM <- suppressMessages(melt(pctM))

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
  geom_point() + suppressMessages(geom_smooth(method="lm")) +
  labs(title = "Diff Methylation vs Log2 Fold Change",
       y = "% DNA methylation increase in RL-VZ",
       x = "Log2 Fold Change in RL-VZ") +
  theme_minimal(base_size = 12)  +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_vline(xintercept = 0, linetype="dashed", color = "red") + 
  ggtitle(sprintf("Promoters of RL DEGs (>14 PCW) (upstream %1.1fkb, downstream %1.1fkb)\n(MinCVG=%i, MinSamp=%i)",
    upstream/1000, downstream/1000, minCvg, minSamp))

p <- p + stat_cor(method="pearson",label.y=25, label.x=1)
p <- p + geom_text_repel(
  data = subset(sig, abs(diffMethy) > 0.05 & abs(logFC) > 0.1),
  aes(label = gene),
  size = 3,
  nudge_x = 0.1,
  nudge_y = 0.1,
  show.legend = FALSE
)
ggsave(
  filename = sprintf("%s/AllPromoters_diffMethy_vs_Log2FC_sig_minCvg%i_up%i_down%i.pdf",
    outDir, minCvg, upstream, downstream),
  plot = p,
  width = 8, height = 6
)

cr <- cor.test(sig$logFC, sig$diffMethy, method="pearson")
cat(sprintf("Correlation between promoter methylation change and gene expression change (All Promoters, DEG genes):\n"))
cat(sprintf("N=%i, Pearson's r = %1.3f, p-value = %1.3e\n", 
  nrow(sig), cr$estimate, cr$p.value))




}, error = function(e) {
  cat("Error in setting options: ", e$message, "\n")
}, finally={
  sink()
})







