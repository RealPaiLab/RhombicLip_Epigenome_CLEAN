# relabel Mannens cells to clean the labels clean VZ-SVZ label
rm(list=ls())

library(Seurat)
library(ggplot2)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(AnnotationHub)

outRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024"
inFile <- sprintf("%s/250612/Mannens_2024_seurat.qs", outRoot)

updatedFile <- sprintf("%s/250612/Mannens_2024_seurat.qs", outRoot)

dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/Relabelled_%s", outRoot, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

message("Using existing Seurat object with updated fragment paths.")
  srat <- qs::qread(updatedFile)

genes2test <- c("ATOH1", "WLS","SOX2", "MKI67","LMX1A","EOMES","RBFOX3",
"PAX6","EGR1","NEUROD1","NEUROD2","TBR1","RELN")
cellLevels <- c("RL-VZ","RL-SVZ","UBC precursors","UBC","GCP","GC","Myeloid", "Endothelial")

DefaultAssay(srat) <- "RNA"
ct <- srat[["RNA"]]$counts

srat$celltype_final <- factor(srat$celltype_final, 
    levels = cellLevels)
p1 <- DotPlot(
  srat, features = genes2test,
  group.by = "celltype_final"
)  
ggsave(
  filename = sprintf("%s/Mannens_DotPlot_beforeRelabel_%s.pdf", outDir, dt),
  plot = p1, width = 15, height = 5
)

lbl <- srat$celltype_final
newlbl <- srat$celltype_final
idx <- which(ct["WLS",]>0 & lbl %in% c("RL-SVZ"))
cat(sprintf("* %i cells have WLS>0 & label RL-SVZ. \nRelabelling to RL-VZ\n", 
    length(idx)))

idx <- which(ct["SOX2",]>0 & lbl %in% c("RL-SVZ"))
cat(sprintf("* %i cells have SOX2>0 & label RL-SVZ. \nRelabelling to RL-VZ\n", 
    length(idx)))
newlbl[idx] <- "RL-VZ"

idx <- which(ct["SOX2",]<0 & ct["WLS",]<0 & lbl %in% c("RL-VZ"))
cat(sprintf("* %i cells have no SOX2 or WLS. \nRelabelling to RL-VZ\n", 
    length(idx)))
newlbl[idx] <- "RL-SVZ"


idx <- which(ct["PAX6",]>0 & lbl %in% c("GC"))
cat(sprintf("* %i cells have PAX6>0 & label GC. \nRelabelling to GCP\n", 
    length(idx)))
newlbl[idx] <- "GCP"
idx <- which(ct["PAX6",]>0 & ct["EOMES",]<1 & lbl %in% c("UBC"))
cat(sprintf("* %i cells have PAX6>0, no EOMES but labelled UBC. \nRelabelling to GCP\n", 
    length(idx)))
newlbl[idx] <- "GCP"

idx <- which(ct["PAX6",]>0 & lbl %in% c("UBC precursors", "UBC"))
cat(sprintf("* %i cells have PAX6>0, no EOMES but labelled UBC. \nRelabelling to GCP\n", 
    length(idx)))
newlbl[idx] <- "GCP"

idx <- which(lbl %in% c("UBC precursors"))
cat("* Mapping all UBC precursors to UBC")
newlbl[idx] <- "UBC"


srat$cleaner_cellLabel <- newlbl
srat$cleaner_cellLabel <- factor(
  srat$cleaner_cellLabel, 
  levels = cellLevels
)
p1 <- DotPlot(
  srat, features = genes2test,
  group.by = "cleaner_cellLabel"
) 
ggsave(
  filename = sprintf("%s/Mannens_DotPlot_afterRelabel_%s.pdf", outDir, dt),
  plot = p1, width = 15, height = 5
)

# plot UMAP using SCT assay
DefaultAssay(srat) <- "RNA"
p1 <- DimPlot(srat, group.by="cleaner_cellLabel")
ggsave(
  filename = sprintf("%s/Mannens_UMAP_afterRelabel_%s.pdf", outDir, dt),
  plot = p1,
  width = 6, height = 6
)

write.table(newlbl, 
  file = sprintf("%s/Mannens_relabelled_cellLabels_%s.txt", outDir, dt),
  quote = FALSE, row.names = TRUE, col.names = FALSE)



