rm(list=ls())

library(Seurat)
inDir <- "/home/rstudio/isilon/src/neurodev-genomics/scRNAseq/Zhong_2023"

cat("Reading Zhong et al.\n")
t0 <- Sys.time()
dat <- Read10X(inDir)
print(Sys.time() - t0)

cat("create Seurat object\n")
srat <- CreateSeuratObject(
  counts = dat, 
  project = "Zhong", 
  min.cells = 3,    # Include features detected in at least 3 cells
  min.features = 200 # Include cells with at least 200 detected genes
)

