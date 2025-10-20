# code for motif enrichment in Mannens.
rm(list=ls())

library(Seurat)
library(ggplot2)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(chromVAR)

outRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024"
inFile <- sprintf("%s/250612/Mannens_2024_seurat.qs", outRoot)

updatedFile <- sprintf("%s/Mannens_Motif/Mannens_2024_WithMotifAndChromVAR_%s.qs", 
  outRoot, "250619")

outDir <- sprintf("%s/Mannens_Motif/%s", outRoot, format(Sys.Date(), "%y%m%d"))
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

# plot dotplot of ChromVAR motif activity
# @param obj Seurat object with chromVAR assay
# @param motifs Character vector of motif IDs to plot or gene names if useGeneNames is TRUE
# @param useGeneNames Logical, if TRUE use gene names instead of motif IDs
# @return ggplot object of the dot plot
plotChromVarDotPlot <- function(obj, motifs, plotIdents=NULL){
  if (is.null(plotIdents)) plotIdents <- levels(Idents(obj))

  p <- DotPlot(obj, features = motifs, assay="chromvar",group.by="celltype_final", idents = plotIdents)

  # change from motif id to motif name
  lv <- levels(p$data$features.plot)
  lb <- ConvertMotifID(obj, id=lv)
  p$data$features.plot <- factor(
      p$data$features.plot,
      levels = lv,
      labels = lb
  )
  p 
}

logFile <- sprintf("%s/Mannens_analyze.log", outDir)
sink(logFile,split=TRUE)

tryCatch({
    if (!file.exists(updatedFile)){
        cat("* Going to create a new Seurat object with motifs and chromVAR\n")
        cat("* Reading Mannens file with annotations & fragment files updated\n")
        srat <- qs::qread(inFile)

        DefaultAssay(srat) <- "peaks" # switch to peaks assay
        Idents(srat) <- "celltype_final"
        srat$celltype_final <- factor(
         srat$celltype_final,
          levels = c("RL-VZ", "RL-SVZ", "UBC precursors", "UBC","GCP","GC","Endothelial", "Myeloid")
        )

        cat("* remove non-standard chromosomes\n")
        # otherwise you will get a subscript out of bounds error when running
        # motif enrichment. The motifs
        # database won't have entries for the non-standard chromosomes.
        # so your motifobject and seurat object will have different indexing 🙀
        # a symptom of this is a warning message for the AddMotifs() call.
         gr <- granges(srat)
         names(gr) <- rownames(srat)
        keep_chr <- paste0("chr", c(1:22, "X", "Y")) # Define standard chromosomes
        gr_filt <- gr[seqnames(gr) %in% keep_chr]
        idxtokeep <- which(rownames(srat) %in% names(gr_filt))
        cat("before = ", length(gr), " after = ", length(gr_filt), "\n")
        srat <- srat[idxtokeep, ]

        cat("* Get motifs from JASPAR2020\n")
        # Adding motif information
        pfm <- getMatrixSet(
          x = JASPAR2020,
          opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
        )

        cat("* Add motifs to Seurat objects - takes ~10-12 min\n")
        t0 <- Sys.time()
        srat <- AddMotifs(
          object = srat,
          genome = BSgenome.Hsapiens.UCSC.hg38,
          pfm = pfm
        )
        print(Sys.time() - t0)

        cat("* Process peaks, perform dimensionality reduction\n")
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

        cat("* Running chromVAR to compute motif deviations\n")
        t0 <- Sys.time()
        srat <- RunChromVAR(srat, genome = BSgenome.Hsapiens.UCSC.hg38)
        print(Sys.time() - t0)

        cat("Saving updated Seurat object with motifs and chromVAR\n")
        qs::qsave(srat, file = updatedFile)

    } else {
        cat("* Using existing Seurat object with motifs and chromVAR\n")
        srat <- qs::qread(updatedFile)
    }

DefaultAssay(srat) <- "chromvar"
Idents(srat) <- "celltype_final"
srat$celltype_final <- factor(
  srat$celltype_final,
  levels = c("RL-VZ", "RL-SVZ", "UBC precursors", "UBC","GCP","GC","Endothelial", "Myeloid")
)
p <- DimPlot(srat)
ggsave(p, file=sprintf("%s/dimplot_chromvar.pdf",outDir))

DefaultAssay(srat) <- "peaks"
darPairs <- list(
  c("RL-SVZ", "Myeloid"),
  c("RL-VZ","RL-SVZ")
)

TFsofinterest <- c("OTX2","EOMES","NEUROD1","NEUROD2",
  "Atoh1","ATOH1(var.2)",
  "BARHL1", "SOX2","Mafb","Crx","NRL","EGR1","TP53"
)
DefaultAssay(srat) <- "peaks"
p <- plotChromVarDotPlot(srat, ConvertMotifID(srat, name=TFsofinterest))
p <- p + 
  ggtitle("ChromVAR for rhombic lip TFs: Mannens 2024") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file=sprintf("%s/dotplot_chromvar_tfs.pdf",outDir),
width=10,height=6)


for (pair in darPairs) {
  DefaultAssay(srat) <- "peaks"
  oDir <- sprintf("%s/%s_vs_%s", outDir, pair[1], pair[2])
  if (!dir.exists(oDir)) {
    dir.create(oDir, recursive = FALSE)
  }

  cat(sprintf("\n* Processing pair: %s vs %s\n", pair[1], pair[2]))
  cat("Reading DARs\n")
  daPeakFile <- sprintf("%s/250612/dar_peaks_%s_vs_%s/da_peaks_%s_vs_%s.csv",  
    outRoot,pair[1],pair[2],pair[1],pair[2])
  
  da_peaks <- read.delim(daPeakFile, sep="\t", header=TRUE)
  notthere <- which(!da_peaks$peak %in% rownames(srat))
  cat("Number of peaks not in Seurat object: ", length(notthere), "\n")
  if (length(notthere) > 0) {
    da_peaks <- da_peaks[-notthere, ]
  }

  cat("Get foreground and background peaks\n")
  # get top differentially accessible peaks
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 
    & da_peaks$pct.1 > 0.2, ])
  open.peaks <- AccessiblePeaks(srat, idents = pair)

   meta.feature <- GetAssayData(srat, 
            assay = "peaks", layer = "meta.features")
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    n = 50000
  )

  t0 <- Sys.time()
  enriched.motifs <- FindMotifs(
    object = srat,
    features = top.da.peak,
    background = peaks.matched
  )
  print(Sys.time() - t0)

  # write the enriched motifs to a file
  motifFile <- sprintf("%s/%s_vs_%s_enriched_motifs.csv", oDir, 
    pair[1], pair[2])
  cat("* Writing enriched motifs to file: ", motifFile, "\n")
  write.csv(enriched.motifs, file = motifFile, row.names = TRUE)

  # plot the top enriched motifs
  top.motifs <- head(enriched.motifs, n=15)

#  DefaultAssay(srat) <- "chromvar"
  p <- plotChromVarDotPlot(srat, top.motifs$motif)
  p <- p + 
    ggtitle(sprintf("Top enriched motifs for %s vs %s", pair[1], pair[2])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  ggsave(
    filename = sprintf("%s/top_motifs_%s_vs_%s.pdf", oDir, 
      pair[1], pair[2]),
    plot = p,
    width = 10, height = 6
  )
}

},error=function(ex){
  print(ex)
},finally={
  sink()
})