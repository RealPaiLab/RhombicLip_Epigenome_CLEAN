# Mannens fetal hindbrain scMultiome set. 
# EDA, updated fragment file paths and annotation, and DAR calling.
rm(list=ls())

library(Seurat)
library(ggplot2)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(AnnotationHub)
#library(JASPAR2020)
#library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(patchwork)

inFile <- "/home/rstudio/isilon/src/neurodev-genomics/multiome/Mannens_2024/from_anders/rl_micro_mannens.qs"

chromSizes <- "/home/rstudio/isilon/src/ucsc-goldenpath/hg38/hg38-20221124T1111.chrom.sizes"
geneDef <- "/home/rstudio/isilon/src/gencode/GRCh38/gencode.v42.basic.annotation.gtf"
fragRoot <- "/home/rstudio/isilon/src/neurodev-genomics/multiome/Mannens_2024/from_anders/fragment_files/"
samps <- list(
  "10X280_1_ABCD_1"="10X280_1_ABCD_1/atac_fragments.tsv.gz",
  "10X346_1_ABCD_1"="10X346_1_ABCD_1/atac_fragments.tsv.gz",
  "10X365_2_ABCDE_2"="10X365_2_ABCDE_2/atac_fragments.tsv.gz",
  "10X406_1_ABCD_2"="10X406_1_ABCD_2/atac_fragments.tsv.gz"
)

outRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outRoot, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

updatedFile <- sprintf("%s/250612/Mannens_2024_seurat.qs", outRoot)
newLabels <- sprintf("%s/Relabelled_250703/Mannens_relabelled_cellLabels_250703.txt", outRoot)

logFile <- sprintf("%s/Mannens_analyze.log", outDir)
sink(logFile,split=TRUE)

tryCatch({
if (!file.exists(updatedFile)) {
  cat("* Reading Mannens file from AWE\n")
  srat <- qs::qread(inFile)

  cat("* Updating fragment paths in Seurat object\n")
  frag <- Fragments(srat@assays$peaks)
  Fragments(srat@assays$peaks) <- NULL
  olddir <- "/hpf/largeprojects/mdtaylor/Datasets_from_papers/sc_multiome/mannens_2023_biorxiv/10X_outs_cerebellum_hindbrain/"
  for (k in 1:length(frag)){
    p <- frag[[k]]@path
    p <- sub(olddir,fragRoot, p)
    frag[[k]] <- UpdatePath(frag[[k]], new.path = p)
  }
  cat("About to assign updated fragments\n")
  Fragments(srat@assays$peaks) <- frag

  cat("\n* Annotating Seurat object with gene information\n")
  ah  <- AnnotationHub()
  ensdb_v98 <- ah[["AH75011"]]
  annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "hg38"
  DefaultAssay(srat) <- "peaks"
  Annotation(srat) <- annotations

  cat("Saving updated Seurat object with fragment paths and gene annotation data\n")
  qs::qsave(srat, file = updatedFile)
  message("Seurat object saved with updated fragment paths.")
} else {
  message("Using existing Seurat object with updated fragment paths.")
  srat <- qs::qread(updatedFile)
}

cat("*** Adding new labels for celltype\n")
newlbls <- read.table(newLabels, 
  header = FALSE, stringsAsFactors = FALSE)

md <- srat[[]]
midx <- match(rownames(md), newlbls$V1)
if (all.equal(newlbls$V1[midx], rownames(md)) != TRUE) {
  stop("Mismatch in rownames of metadata and new labels.")
}
srat$cleaner_cellLabel <- newlbls$V2[midx]

cat("HENCEFORTH USING cleaner_cellLabel instead of celltype_final\n")
srat$celltype_final <- srat$cleaner_cellLabel

# print num cells and features
message(sprintf("Number of cells: %d", ncol(srat)))
message(sprintf("Number of features: %d", nrow(srat)))

x <- srat[[]][,c("Donor","sample","SEX","Age")]
x <- x[!duplicated(x),]
message(sprintf("Number of unique donors: %d", nrow(x)))
message(sprintf("Number of unique samples: %d", length(unique(srat$sample))))
print(table(x$SEX,useNA="always"))

cat("Cells by donor\n")
print(table(srat$Donor, useNA="always"))

# process the ATAC peaks
DefaultAssay(srat) <- "peaks"
srat <- FindTopFeatures(srat, min.cutoff = 10)
t0 <- Sys.time()
srat <- RunTFIDF(srat)
t1 <- Sys.time()
message(sprintf("Time taken for TF-IDF: %s", t1 - t0))
t0 <- Sys.time()
srat <- RunSVD(srat)
t1 <- Sys.time()
message(sprintf("Time taken for SVD: %s", t1 - t0))
srat <- RunUMAP(srat, reduction = "lsi", dims=2:30)

p1 <- DimPlot(srat, group.by="celltype_final")
ggsave(
  filename = sprintf("%s/umap_peaks.pdf", outDir),
  plot = p1,
  width = 6, height = 6
)

DefaultAssay(srat) <- "SCT"
p1 <- DimPlot(srat, group.by="celltype_final")
ggsave(
  filename = sprintf("%s/umap_sct.pdf", outDir),
  plot = p1,
  width = 6, height = 6
)

DefaultAssay(srat) <- "RNA"
p1 <- DotPlot(
  srat, features = c("ATOH1", "WLS","MKI67","EOMES","LMX1A","RBFOX3"),
    #"OTX2","PTPRZ1","HOPX","DCX","NEUROD1","ASCL1","TUBB3"),
  group.by = "celltype_final"
)
ggsave(
  filename = sprintf("%s/dotplot_celltype_markers.pdf", outDir),
  plot = p1,
  width = 16, height = 6
)

p1 <- DimPlot(srat, group.by="celltype_final")
ggsave(
  filename = sprintf("%s/umap_rna.pdf", outDir),
  plot = p1,
  width = 6, height = 6
)

plotCvg <- function(rg, up=20, down=20,pfx="",ttl="",oDir){
  x <- unlist(strsplit(rg, "-"))
  gr <- GRanges(x[1],IRanges(as.numeric(x[2]), as.numeric(x[3])))
  p <- CoveragePlot(
    object = srat,
    region = rg,
    extend.upstream = up*1000,
    extend.downstream = down*1000,
    region.highlight = gr
  )
  p <- p + ggtitle(ttl)
  ggsave(
    filename = sprintf("%s/%s_%s_coverage.pdf", oDir, pfx, rg),
    plot = p,
    width = 6, height = 4
  )
}

# call DARs and annotate for two groups.
getAnnotateDARs <- function(g1,g2){
  message(sprintf("Finding DA peaks between %s and %s", g1, g2))

  darDir <- sprintf("%s/dar_peaks_%s_vs_%s", outDir, g1, g2)
  if (!dir.exists(darDir)) {
    dir.create(darDir, recursive = FALSE)
  }

cat("Calling FindMarkers with test.use='LR', min.pct=0.05, and latent.vars='nCount_peaks'\n")
  t0 <- Sys.time()
  da_peaks <- FindMarkers(
    object = srat, 
    ident.1 = g1,
    ident.2 = g2,
    test.use="LR",
    min.pct = 0.05,
    latent.vars = "nCount_peaks"
  )

  da_peaks$peak <- rownames(da_peaks)
  da_peaks$celltype1 <- g1
  da_peaks$celltype2 <- g2
  write.table(da_peaks, 
    file = sprintf("%s/da_peaks_%s_vs_%s.csv", darDir, g1, g2),
    sep="\t",col=TRUE,row=TRUE,
    quote = FALSE
  )
  print(sprintf("Time taken to find DA peaks: %s", Sys.time() - t0))
  cat("\n Num significant DA peaks: ", sum(da_peaks$p_val_adj < 0.05))

  # plot QQ plot
  pdf(sprintf("%s/da_peaks_qqplot_%s_vs_%s.pdf", darDir, g1, g2), width = 6, height = 6)
  hist(da_peaks$p_val, breaks=200, 
       main = sprintf("Histogram of p-values for DA peaks: %s vs %s", g1, g2),
       xlab = "p-value", ylab = "Frequency")
  dev.off()

  # plot volcano plot
  pdf(sprintf("%s/volcano_%s_vs_%s.pdf", darDir, g1, g2), width = 6, height = 6)
  ggplot(da_peaks, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.3), alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = sprintf("Volcano plot of DA peaks: %s vs %s", g1, g2),
         x = "Log2 Fold Change", y = "-log10 Adjusted p-value")
  dev.off()
  
  # plot coverage for top peaks
  tmp <- subset(da_peaks, p_val_adj < 0.05 & avg_log2FC > 2)
  nr <- min(5,nrow(tmp))
  for (k in 1:nr){
    ttl <- sprintf("%s vs %s: %s", g1, g2, rownames(tmp)[k])
    plotCvg(rownames(tmp)[k],
      sprintf("%s-%s",g1,g2),ttl=ttl,up=40,down=40,oDir=darDir)
  }
}

DefaultAssay(srat) <- "peaks"
Idents(srat) <- "celltype_final"
#getAnnotateDARs("RL-VZ", "RL-SVZ")
#getAnnotateDARs("RL-VZ","UBC")
#getAnnotateDARs("RL-VZ","GCP")
#getAnnotateDARs("RL-SVZ", "Myeloid")
#getAnnotateDARs("RL-SVZ","UBC precursors")
#getAnnotateDARs("UBC precursors","UBC")
#getAnnotateDARs("RL-SVZ","GCP")
#getAnnotateDARs("GCP","GC")

# link peaks to genes
cat("\n* Linking peaks to genes\n")
srat <- RegionStats(srat, genome = BSgenome.Hsapiens.UCSC.hg38)
t0 <- Sys.time()
srat <- LinkPeaks(
  object = srat,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = rownames(srat[["SCT"]]),
  min.cells = 100
  #genes.use = "SOX2"
)
print(sprintf("Time taken to link peaks to genes: %s", Sys.time() - t0))

linkGR <- Links(srat)
# save the linked peaks
write.table(as.data.frame(linkGR),
  file = sprintf("%s/Mannens_linked_peaks_%s.txt", outDir,dt),
  sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE
)

p1 <- CoveragePlot(
  srat, region = "SOX2", features = "SOX2", expression.assay = "SCT",
  extend.upstream = 2*1000, extend.downstream = 10*1000
)
ggsave(
  filename = sprintf("%s/coverage_SOX2.pdf", outDir),
  plot = p1,
  width = 8, height = 4
)

browser()




#### plot peaks around genes of interest
###gene_anno <- rtracklayer::import(geneDef)
###gene_df <- as.data.frame(gene_anno)
###gene_df <- subset(gene_df, type=="gene" & gene_type=="protein_coding")
###getTSS <- function(nm, up=10,down=10){
###  x <- subset( gene_df, gene_name == nm )
###  if (x$strand == "+")
###    rg <- paste(x$seqnames,x$start, x$start,sep="-")
###  else 
###    rg <- paste(x$seqnames,x$end, x$end,sep="-")
###rg 
###  p <- CoveragePlot(
###    object=srat, region=rg,
###    extend.upstream=up*1000, 
###    extend.downstream=down*1000
###  )
###  ggsave(p, 
###    file=sprintf("%s/%s_Coverage.pdf",outDir,nm),
###    width=10,height=6)
###}
###
###for (g in c("MKI67","RBFOX3")){#c("WLS","SOX2","NEUROD1","BRCA1",
###  print(g)
###  getTSS(g)
###}


}, error = function(e) {
  message("Please install the 'qs' package to read the input file.")
}, finally= {
  sink()
}
)
