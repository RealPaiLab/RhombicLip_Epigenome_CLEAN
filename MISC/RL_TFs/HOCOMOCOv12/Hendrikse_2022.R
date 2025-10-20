# Used this script to obtain a list of gene that is likely expressed in RL cells.
# To be robust, used a very low threshold - i.e. expressed in at least 1% cells in one cell type.
# Got gene list for VZ and SVZ separately to avoid extreme condition - gene expressed in 1% VZ cells 
# but 0 in SVZ; if pool together the gene will be lost. After obtaining the two gene lists, they were
# combined to give genes expressed in RL. This list can be used to remove TFs not even expressed in
# the RL cells to reduce multiple testing burden and interpretation complexity. 

library(Seurat)

file_dir <- "/data/xsun/input/Hendrikse_2022/glutamatergic_dev_Liam.RDS"
out_dir <- "/data/xsun/db/meme/HOCOMOCOv12/activeGenes"

### Load data
## snRNA-seq
message("Loading Hendrikse 2022 glutamatergic lineage snRNA-seq dataset")
glu <- readRDS(file_dir)
glu@active.assay <- "RNA"


### Subset RL-VZ/SVZ populations
message("Subsetting RL-VZ cells")
vz_cells <- subset(glu, idents = "RL-VZ")
print(vz_cells)

message("Subsetting RL-SVZ cells")
svz_cells <- subset(glu, idents = "RL-SVZ")
print(svz_cells)


### Check expression
thresh <- 0.01 # set gene should be expressed in at least 1% of cells

## RL-VZ
vz_counts <- GetAssayData(object = vz_cells, slot = "counts")
vz_nonzero <- rowSums(vz_counts > 0)
vz_prop <- vz_nonzero/ncol(vz_counts)
vz_genes <- names(vz_prop)[vz_prop >= thresh]
message(sprintf("%s genes are expressed in at least %.2f%% cells in RL-VZ", length(vz_genes), thresh*100))

## RL-SVZ
svz_counts <- GetAssayData(object = svz_cells, slot = "counts")
svz_nonzero <- rowSums(svz_counts > 0)
svz_prop <- svz_nonzero/ncol(svz_counts)
svz_genes <- names(svz_prop)[svz_prop >= thresh]
message(sprintf("%s genes are expressed in at least %.2f%% cells in RL-SVZ", length(svz_genes), thresh*100))

### output list
rl_genes <- unique(c(vz_genes, svz_genes))
message(sprintf("%s genes are expressed in at least %.2f%% cells in RL-VZ or RL-SVZ", length(rl_genes), thresh*100))

vz_outFile <- sprintf("%s/Hendrikse2022_RLVZ_activeGenes", out_dir)
svz_outFile <- sprintf("%s/Hendrikse2022_RLSVZ_activeGenes", out_dir)
rl_outFile <- sprintf("%s/Hendrikse2022_RL_activeGenes", out_dir)

write.table(vz_genes, vz_outFile, col.names = F, row.names = F, quote = F)
write.table(svz_genes, svz_outFile, col.names = F, row.names = F, quote = F)
write.table(rl_genes, rl_outFile, col.names = F, row.names = F, quote = F)




