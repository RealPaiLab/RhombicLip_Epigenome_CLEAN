rm(list=ls())

# Infer enhancer-gene peaks from Sarropoulos ATAC-seq data
# https://stuartlab.org/signac/articles/pbmc_multiomic

library(Seurat)
library(SeuratWrappers)
library(Signac)
library(SingleCellExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationHub)
library(ggplot2)
###library(JASPAR2020)
###library(TFBSTools)
library(cicero)

inDir <- "/home/rstudio/isilon/src/neurodev-genomics/scMultiome/Sarropoulos_2026/"
outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/"

processedDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260402/"
atacRDS <- paste0(processedDir, "Sarropoulos_ATAC_Seurat.rds")
rnaRDS <- paste0(processedDir, "Sarropoulos_RNA_Seurat.rds")

geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if(!dir.exists(outDir)) dir.create(outDir, recursive=FALSE)

logFile <- sprintf("%s/%s_Sarropoulos_ATAC.log", outDir, dt)
sink(logFile, append=TRUE)

tryCatch({

if (file.exists(atacRDS) & file.exists(rnaRDS)){
  cat("Loading existing Seurat objects...\n")
  srat <- readRDS(atacRDS)
  rna <- readRDS(rnaRDS)
} else {
  cat("Processing Sarropoulos ATAC and RNA data from scratch...\n")
  atacFile <- paste0(inDir, "human_atac_sce.rds")
  rnaFile <- paste0(inDir, "human_rna_sce.rds")
  #mdataFile <- paste0(inDir, "science.adw9154_tables_s1_to_s18.xlsx")

  cat("reading Sarropoulos ATAC file\n")
  t0 <- Sys.time()
  atac <- readRDS(atacFile)
  print(Sys.time()-t0)
  cat(sprintf("Sarropoulos ATAC: %i cells, %i features\n", ncol(atac), nrow(atac)))

  counts <- assay(atac)
  rng <- sub("hg38_","", rownames(counts))
  cpos <- regexpr(":", rng); 
  hpos <- regexpr("-", rng);
  chrom <- substr(rng, 1, cpos-1)
  start <- as.integer(substr(rng, cpos+1, hpos-1))
  end <- as.integer(substr(rng, hpos+1, nchar(rng)))
  peaks <- GRanges(seqnames=chrom, ranges=IRanges(start=start, end=end))

  mdata <- colData(atac)
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    ranges = peaks,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = NULL
  )

  srat <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = mdata
  )

  cat("Only keep mid-gestation samples to reduce compute time\n")
  cells.tokeep <- which(mdata$Stage_group %in% c("11wpc","15-17wpc"))
  cat(sprintf("Keeping %i cells\n", length(cells.tokeep)))
  srat <- srat[, cells.tokeep]
  mdata <- mdata[cells.tokeep, ]

  if (all.equal(rownames(mdata), colnames(srat))!=TRUE){
    cat("fix the ordering of metadata and srat\n")
  }

  srat$rna_precisest_label <- mdata$rna_precisest_label
  srat$rna_dev_stage <- mdata$rna_dev_state
  srat$Stage_group <- mdata$Stage_group

  ###### all peaks are from standard chroms chr1-X.
  ######peaks.keep <- seqnames(granges(srat)) %in% standardChromosomes(granges(srat))
  ######srat <- srat[as.vector(peaks.keep), ]
  #####
  ###### QC -- the commented-out lines could not be run because the data were not provided with 
  ###### a fragment file.
  ###### compute nucleosome signal score per cell
  ######srat <- NucleosomeSignal(object = srat)
  ###### compute TSS enrichment score per cell
  ######srat <- TSSEnrichment(object = srat)
  ###### add fraction of reads in peaks
  ######srat$pct_reads_in_peaks <- srat$peak_region_fragments / srat$passed_filters * 100
  #####
  ###### add blacklist ratio
  #####blacklist_regions <- ah[['AH107305']] # blacklist regions for hg38
  #####srat$blacklist_ratio <- FractionCountsInRegion(
  #####  object = srat, 
  #####  assay = 'peaks',
  #####  regions = blacklist_regions
  #####)
  #####
  #####cat("Blacklist ratio is very low, indicating good data quality\n")
  #####print(summary(srat$blacklist_ratio))
  ########> summary(srat$blacklist_ratio)
  ########     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
  ########0.000e+00 0.000e+00 0.000e+00 4.317e-05 5.476e-05 1.999e-03
  #####
  cat("Running dimensionality reduction...\n")
  srat <- RunTFIDF(srat)
  srat <- FindTopFeatures(srat, min.cutoff = 'q0')
  t0 <- Sys.time()
  srat <- RunSVD(srat)
  print(Sys.time()-t0)

  cat("Finding clusters...\n")
  srat <- RunUMAP(object = srat, reduction = 'lsi', dims = 2:30)
  srat <- FindNeighbors(object = srat, reduction = 'lsi', dims = 2:30)
  srat <- FindClusters(object = srat, verbose = FALSE, algorithm = 3)
  p <- DimPlot(srat, label=TRUE) + ggtitle("Sarropoulos ATAC clusters")
  ggsave(sprintf("%s/DimPlot.png", outDir), plot=p, width=6, height=5)

  p <- DimPlot(srat, group.by="rna_precisest_label", label=TRUE) + ggtitle("Sarropoulos ATAC clusters colored by RNA labels")
  ggsave(sprintf("%s/DimPlot_by_rna_precisest_label.png", outDir), plot=p, width=14, height=8)

  p <- DimPlot(srat, group.by="rna_precisest_label", label = FALSE) + ggtitle("Sarropoulos ATAC clusters colored by RNA labels")
  ggsave(sprintf("%s/DimPlot_by_rna_precisest_label_nolab.png", outDir), plot=p, width=14, height=8)

  # DimPlot showing only the progenitor_RL cluster
  p <- DimPlot(srat, group.by="rna_precisest_label", label = FALSE)
  p <- p + ggtitle("Sarropoulos ATAC clusters colored by RNA labels") +
    theme(legend.text = element_text(size=6)) +
    scale_color_manual(values=c("progenitor_RL"="red", "other"="grey"))
  ggsave(sprintf("%s/DimPlot_by_rna_precisest_label_progenitor_RL.png", outDir), plot=p, width=14, height=8)

saveRDS(srat, file=sprintf("%s/Sarropoulos_ATAC_Seurat.rds", outDir))
  
  cat("------------------------------------\n")
  cat("reading Sarropoulos RNA file\n")
  t0 <- Sys.time()
  rna <- readRDS(rnaFile)
  print(Sys.time()-t0)
  cat(sprintf("Sarropoulos RNA: %i cells, %i features\n", ncol(rna), nrow(rna)))
  rna_sce <- as(rna, "SingleCellExperiment")
  ###
  rna <- as.Seurat(rna_sce, counts="counts", data = NULL)
  rna <- subset(rna, subset = Stage_group %in% c("11wpc","15-17wpc"))
  cat(sprintf("Keeping %i cells\n", ncol(rna)))

  options(future.globals.maxSize= 100 * 1024^3)
  t0 <- Sys.time()
  rna <- SCTransform(rna, assay="originalexp",
    ncells = 5000, conserve.memory = TRUE, 
    return.only.var.genes = FALSE)
  print(Sys.time()-t0)

  rna <- RunPCA(rna, assay="SCT", npcs=50)
  eb <- ElbowPlot(rna, ndims=50)
  ggsave(sprintf("%s/RNA_ElbowPlot.png", outDir), plot=eb, width=6, height=5, bg="white")
  rna <- RunUMAP(rna, dims=1:30)

  p <- DimPlot(rna, group.by="precisest_label", label=FALSE) + ggtitle("Sarropoulos RNA clusters")
  ggsave(sprintf("%s/RNA_DimPlot_by_precisest_label.png", outDir), plot=p, width=18, height=8)

  # Show a DimPlot only highlighting progenitor_RL and progenitor_RL_early
  p <- DimPlot(rna, group.by="precisest_label", label=FALSE)
  p <- p + ggtitle("Sarropoulos RNA clusters colored by RNA labels") +
    theme(legend.text = element_text(size=6)) +
    scale_color_manual(values=c("progenitor_RL"="red", "progenitor_RL_early"="orange", "other"="grey"))
  ggsave(sprintf("%s/RNA_DimPlot_by_precisest_label_progenitor_RL.png", outDir), plot=p, width=6, height=6)

  saveRDS(rna, file=sprintf("%s/Sarropoulos_RNA_Seurat.rds", outDir))
}

cat("Adding gene annotations for GeneActivity below.\n")
ah <- AnnotationHub()
# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(srat) <- annotations
#####
###### plot a proportion barplot of rna_precisest_label within each ATAC cluster
#####md <- srat[[]]
#####p <- ggplot(md, aes(x=seurat_clusters, fill=rna_precisest_label)) +
#####  geom_bar(position="fill") +
#####  xlab("ATAC cluster") +
#####  ylab("Proportion of cells") +
#####  ggtitle("Proportion of RNA labels within each ATAC cluster") +
#####  theme(legend.text = element_text(size=6))
#####ggsave(sprintf("%s/rna_precisest_label_by_ATAC_cluster.png", outDir), plot=p, width=14, height=8) 
#####
#####Idents(srat) <- srat$seurat_clusters
#####cat("* Finding marker peaks for the RL progenitor cluster...\n")
#####da_peaks <- FindMarkers(
#####  object = srat,
#####  ident.1 = "20",
#####  ident.2 = NULL,
#####  min.pct = 0.1,
#####  test.use = 'wilcox'
#####)
#####
#####
###cat("* Now adding motif annotation\n")
###pfm <- getMatrixSet(
###  x = JASPAR2020,
###  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
###)
###t0 <- Sys.time()
###srat <- AddMotifs(
###  object = srat,
###  genome = BSgenome.Hsapiens.UCSC.hg38,
###  pfm = pfm
###)
###print(Sys.time()-t0)
###
###srat <- RunChromVAR(
###  object = srat,
###  genome = BSgenome.Hsapiens.UCSC.hg38
###)
###DefaultAssay(srat) <- "chromvar"
###tfs <- c("SOX2","SKOR2","RBFOX2","ZIC1")
###for (g in tfs){
###  p <- FeaturePlot(srat, 
###    features = g,
###    min.cutoff = 'q10',
###    max.cutoff = 'q90',
###    pt.size = 0.1)
###  ggsave(p, filename = sprintf("%s/ChromVAR_%s.png", outDir, g), width=6, height=5)
###}
###
###

# Let's use Cicero to find co-accessible peaks and link them to genes. 
srat.cds <- as.cell_data_set(srat)
srat.cicero <- make_cicero_cds(srat.cds, 
  reduced_coordinates = reducedDims(srat.cds)$UMAP)
srat.cicero <- run_cic

### Cannot integrate snRNAseq-snATACseq as fragment files were not provided for the snATACseq data
#####activity <- GeneActivity(
#####  object = srat,
#####  assay = "peaks",
#####  extend.upstream = 2000,
#####  extend.downstream = 0
#####)
#####
#####
#####rna <- FindVariableFeatures(rna, nfeatures = 5000)
#####transfer.anchors <- FindTransferAnchors(
#####  reference = rna, 
#####  query = srat, 
#####  reduction = 'cca', 
#####  normalization.method = "SCT", 
#####  dims = 1:30
#####)

}, error=function(ex){
  print(ex)
}, finally={
  sink(NULL)
})