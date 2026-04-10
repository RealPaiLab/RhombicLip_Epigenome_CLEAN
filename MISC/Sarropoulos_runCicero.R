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
outDir <- sprintf("%s/%s_Cicero", outDir, format(Sys.Date(),"%y%m%d"))

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

genome <- seqlengths(srat)
genome.df <- data.frame("chr" = names(genome), "length" = genome)

tmp <- rownames(srat)
chroms <- substr(tmp, 1, regexpr("-", tmp)-1)
lv <- levels(factor(chroms))
for (chr in lv) {
    cat(sprintf("Processing %s...\n", chr))
    t0 <- Sys.time()

    # subset seurat object to chromosome and make cicero cds
    srat_chr <- srat[which(chroms == chr),]
    cat(sprintf("Subsetted to %i cells and %i features\n", ncol(srat_chr), nrow(srat_chr)))

    cds <- as.cell_data_set(srat_chr)
    chr.cicero <- make_cicero_cds(cds,
        reduced_coordinates = reducedDims(cds)$UMAP
    )
    print(Sys.time()-t0)
    cat("Running Cicero for chromosome ", chr, "...\n")
    t0 <- Sys.time()
    cur <- run_cicero(chr.cicero, genomic_coords = genome.df[chr, , drop=FALSE])
    print(Sys.time()-t0)

 #   cur <- ciceroList[[chr]]
    #write.table(cur, file=sprintf("%s/cicero_%s.txt", outDir, chr), sep="\t", quote=FALSE, row.names=FALSE)
    save(cur, file=sprintf("%s/cicero_%s.Rdata", outDir, chr))
}

cat("Running Cicero...\n")
t0 <- Sys.time()
conns <- run_cicero(srat.cicero)
print(Sys.time()-t0)

}, error=function(ex){
    print(ex)
}, finally={
    sink(NULL)
})

