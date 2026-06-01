rm(list=ls())

library(Seurat)

rnaFile_wo <-  "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260402/Sarropoulos_RNA_Seurat.rds"
rnaFile_cc <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260416_checkPeaks/Sarropoulos_RNA_Seurat_cellCycleScored.rds"


cat("loading SCT without cell cycle regression")
rna_wo <- readRDS(rnaFile_wo)
DefaultAssay(rna_wo) <- "SCT"
Idents(rna_wo) <- rna_wo$precisest_label
cat("running deg analysis")
deg <- FindMarkers(
    object = rna_wo,
    ident.1 = "progenitor_RL",
    ident.2 = NULL,
    min.pct = 0.1,
    test.use = 'wilcox'
)

cat("loading SCT with cell cycle regression")
rna_cc <- readRDS(rnaFile_cc)
DefaultAssay(rna_cc) <- "SCT"
Idents(rna_cc) <- rna_cc$precisest_label
cat("running deg analysis")
deg_cc <- FindMarkers(
    object = rna_cc,
    ident.1 = "progenitor_RL",
    ident.2 = NULL,
    min.pct = 0.1,
    test.use = 'wilcox'
)

upreg_deg <- deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.1, ]
upreg_deg_cc <- deg_cc[deg_cc$avg_log2FC > 0 & deg_cc$p_val_adj < 0.1, ]

jac <- length(intersect(rownames(upreg_deg), rownames(upreg_deg_cc)))/length(union(rownames(upreg_deg), rownames(upreg_deg_cc)))

