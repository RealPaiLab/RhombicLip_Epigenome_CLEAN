library(Seurat)

# this rds was all MB subtypes scRNA merged tumours using fastMNN
merged_Valdoiu2019_file <- "/Users/Xinghan/data/xsun/output/MISC/Vladoiu2019_merged.rds"
out_dir <- "./output"

### Load file
message("Loading Vladoiu 2019 merged MB scRNA-seq data")
merged_MB <- readRDS(merged_Valdoiu2019_file)
print(merged_MB)

### Drop likely non-tumour clusters (7, 9)
message("Removing likely non-tumour cells")
# DimPlot(merged_MB, label = T, repel = T)
filtered_MB <- subset(merged_MB, idents = c(seq(0,6), 8))
print(filtered_MB)
message(sprintf("%d cells removed", ncol(merged_MB) - ncol(filtered_MB)))

### Subset G3/G4 populations
print(unique(filtered_MB$sample))
## G3
message("Subsetting G3 MB")
G3_MB <- subset(filtered_MB, group == "G3")
print(G3_MB)

G3_samples <- unique(G3_MB$sample)

## G4
message("Subsetting G4 MB")
G4_MB <- subset(filtered_MB, group == "G4")
print(G4_MB)
G4_samples <- unique(G4_MB$sample)


### Obtain active genes
thresh <- 0.01 # set gene should be expressed in at least 1% of cells
## G3
G3_genes <- NULL
for (i in G3_samples) {
  sample_cells <- subset(G3_MB, sample == i)
  sample_counts <- GetAssayData(object = sample_cells, slot = "counts")
  sample_nonzero <- rowSums(sample_counts > 0)
  sample_prop <- sample_nonzero/ncol(sample_counts)
  sample_genes <- names(sample_prop)[sample_prop >= thresh]
  message(sprintf("%s genes are expressed in at least %.2f%% cells in %s", length(sample_genes), thresh*100, i))
  
  G3_genes <- append(G3_genes, sample_genes)
}
G3_genes <- unique(G3_genes)
message(sprintf("%s active genes combined in G3 MB", length(G3_genes)))

## G4
G4_genes <- NULL
for (i in G4_samples) {
  sample_cells <- subset(G4_MB, sample == i)
  print(sample_cells)
  sample_counts <- GetAssayData(object = sample_cells, slot = "counts")
  sample_nonzero <- rowSums(sample_counts > 0)
  sample_prop <- sample_nonzero/ncol(sample_counts)
  sample_genes <- names(sample_prop)[sample_prop >= thresh]
  message(sprintf("%s genes are expressed in at least %.2f%% cells in %s", length(sample_genes), thresh*100, i))
  
  G4_genes <- append(G4_genes, sample_genes)
}

G4_genes <- unique(G4_genes)
message(sprintf("%s active genes combined in G4 MB", length(G4_genes)))


### output list
G34_genes <- unique(c(G3_genes, G4_genes))
message(sprintf("%s active genes combined in G34 MB", length(G34_genes)))

g3_outFile <- sprintf("%s/Vladoiu2019_G3_activeGenes", out_dir)
g4_outFile <- sprintf("%s/Vladoiu2019_G4_activeGenes", out_dir)
g34_outFile <- sprintf("%s/Vladoiu2019_G34_activeGenes", out_dir)

write.table(G3_genes, g3_outFile, col.names = F, row.names = F, quote = F)
write.table(G4_genes, g4_outFile, col.names = F, row.names = F, quote = F)
write.table(G34_genes, g34_outFile, col.names = F, row.names = F, quote = F)





