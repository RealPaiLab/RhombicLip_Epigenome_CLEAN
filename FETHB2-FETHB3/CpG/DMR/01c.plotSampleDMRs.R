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

rnaDEGfile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/251210/edgeR_RLVZvsSVZ_OlderThan14PCW_251210.txt"

olFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs_cCRE_overlap_240712.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_plotRegion"

dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

upstream <- 4000
downstream <- 1000

promFile <- sprintf("%s/promoter_methylation_up%i_down%i.RData", 
  outDir, upstream, downstream)

logFile <- sprintf("%s/plotSampleDMRs.log", outDir)
sink(logFile, split=TRUE)
tryCatch({
pheno <- read.delim(phenoFile, 
        header = TRUE, 
        stringsAsFactors = FALSE
)
dpos <- regexpr("-",pheno$library_ID)
pheno$sample <- substr(pheno$library_ID, 1, dpos-1)


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

# define colours for plots
###if (all.equal(pheno$sample, sampleNames(bsObj$bs))!=TRUE) {
###  stop("Sample names in pheno and bsObj do not match!")
###}
###cols <- c("#7B3294","#008837")[factor(pheno$ROI,
###  levels=c("VZ","SVZ"))]

#### plot detailed methylation view for selected genes
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
    outFile = sprintf("%s/%s_up%i_down%i.pdf", outDir,g, ups[i], downs[i])
  )
}

}, error = function(e) {
  cat("Error in setting options: ", e$message, "\n")
})







