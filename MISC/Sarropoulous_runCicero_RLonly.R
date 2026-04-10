# Run Cicero to infer co-accessible peaks and link them to genes.
rm(list=ls())
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)

atacFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260402/Sarropoulos_ATAC_Seurat.rds"
outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026"
outDir <- sprintf("%s/%s_CiceroRLonly", outDir, format(Sys.Date(),"%y%m%d"))

if (!file.exists(outDir)) dir.create(outDir, recursive = FALSE)
logFile <- sprintf("%s/%s_Sarropoulos_Cicero.log", outDir, format(Sys.Date(),"%y%m%d"))
sink(logFile, split = TRUE)

tryCatch({

cat("reading Sarropoulos ATAC Seurat object\n")
t0 <- Sys.time()
srat <- readRDS(atacFile)
print(Sys.time()-t0)
cat(sprintf("Sarropoulos ATAC: %i cells, %i features\n", ncol(srat), nrow(srat)))

p <- DimPlot(srat, group.by="rna_precisest_label")
ggsave(p, file=sprintf("%s/Sarropoulos_ATAC_UMAP.png", outDir), width=10, height=5)

srat <- subset(srat, subset = rna_precisest_label %in%  c("progenitor_RL"))
cat(sprintf("Subsetted to %i cells and %i features\n", ncol(srat), nrow(srat)))

cds <- as.cell_data_set(srat)
cat("Making Cicero CDS...\n")
rl.cicero <- make_cicero_cds(cds,
    reduced_coordinates = reducedDims(cds)$UMAP
)
genome <- seqlengths(srat)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
genome.df <- subset(genome.df, chr %in% paste0("chr", c(1:22, "X", "Y")))

t0 <- Sys.time()
cicero.res <- run_cicero(rl.cicero, genomic_coords = genome.df)
print(Sys.time()-t0)
save(cicero.res, file=sprintf("%s/cicero_RLonly.Rdata", outDir))

cicero.good <- subset(cicero.res, coaccess > 0.25)
cat(sprintf("Total %i connections with coaccess > 0.25\n", nrow(cicero.good)))
save(cicero.good, file=sprintf("%s/ciceroRLOnly_passCutoff.Rdata", outDir))

}, error=function(ex){
    print(ex)
}, finally={
    sink(NULL)
})

