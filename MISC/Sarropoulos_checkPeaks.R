rm(list=ls())

library(Seurat)
library(Signac)
library(ggplot2)

atacFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260402/Sarropoulos_ATAC_Seurat.rds"
rnaFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260402/Sarropoulos_RNA_Seurat.rds"
outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s_checkPeaks", outDir, dt)

if (!dir.exists(outDir)) dir.create(outDir, recursive=FALSE)

logFile <- sprintf("%s/%s_Sarropoulos_ATAC_checkPeaks.log", outDir, dt)
sink(logFile, split=TRUE)


tryCatch({
cat("reading file...\n")
t0 <- Sys.time()
atac <- readRDS(atacFile)
t1 <- Sys.time()
cat("file read in", t1 - t0, "seconds\n")

pairs <- list(
    c("progenitor_RL", "GCP/UBCP"),
    c("GCP/UBCP", "UBC_Hcrtr2"),
    c("GCP/UBCP", "UBC_Trpc3"),
    c("GCP/UBCP", "GC_diff_2_early")
)
lv <- unlist(lapply(pairs, function(x) paste(x, collapse="_vs_")))

Idents(atac) <- atac$rna_precisest_label
cat("* Finding marker peaks for the RL progenitor cluster...\n")
res <- list()
for (i in seq_along(pairs)){
    cat(sprintf("* Finding marker peaks for %s vs %s...\n", pairs[[i]][1], pairs[[i]][2]))
da_peaks <- FindMarkers(
  object = atac,
  ident.1 = pairs[[i]][1],
  ident.2 = pairs[[i]][2],
  min.pct = 0.1,
  test.use = 'wilcox'
)
res[[i]] <- da_peaks
}

cat("-----------------\n")
cat("ATAC data\n")
cat("-----------------\n")
dirChange <- list()
log2FC <- list()
cat("--------\n")
for (i in seq_along(res)){
    sig <- res[[i]][res[[i]]$p_val_adj < 0.05, ]
    cat(sprintf("%s vs %s\n", pairs[[i]][1], pairs[[i]][2]))
    cat(sprintf("Num significant peaks: %d\n", nrow(sig)))

    print(summary(sig$avg_log2FC))
    cat(sprintf("Num open peaks in group 1: %d\n", sum(sig$avg_log2FC > 0)))
    cat(sprintf("Num open peaks in group 2: %d\n", sum(sig$avg_log2FC < 0)))
    dirChange[[i]] <- cbind(sum(sig$avg_log2FC > 0), sum(sig$avg_log2FC < 0))
    nm <- paste(pairs[[i]], collapse="_vs_")
    log2FC[[nm]] <- sig$avg_log2FC
    cat("--------\n")

}
dirChange <- do.call(rbind, dirChange)
rownames(dirChange) <- sapply(pairs, function(x) paste(x, collapse=" vs "))
colnames(dirChange) <- c("open in group 1", "open in group 2")
print(dirChange)

log2FC_df <- lapply(names(log2FC), function(x) data.frame(pair=x, log2FC=log2FC[[x]]))
log2FC_df <- do.call(rbind, log2FC_df)
log2FC_df$pair <- factor(log2FC_df$pair, levels=lv)
# show violin plots of log2 fold changes for each pair
# add a colour scheme that's sequential for the pairs, e.g. using RColorBrewer
p <- ggplot(log2FC_df, aes(x=pair, y=log2FC, fill=pair)) +
    geom_violin() + 
    geom_boxplot(width=0.1, outlier.shape=NA) +
    scale_fill_brewer(palette="Blues") +
    xlab("log2 fold change") +
    ggtitle("Distribution of log2 fold changes for significant peaks") +
    theme_bw() +
    theme(axis.text=element_text(size=24), legend.position = "none")
ggsave(sprintf("%s/log2FC_density.pdf", outDir), plot=p, width=8, height=5)

# run WmW test for significantly lt 0 peaks in lv[1] and sig gt 0 for lv[3] and lv[4]
t1 <- t.test(log2FC[[lv[1]]], mu = 0, alternative="less")
t3 <- t.test(log2FC[[lv[3]]], mu = 0, alternative="greater")
t4 <- t.test(log2FC[[lv[4]]], mu = 0, alternative="greater")
cat("Significance tests for overall directionality of peak changes\n")
cat(sprintf("%s: mean log2FC = %.3f, p-value = %.3e\n", lv[1], mean(log2FC[[lv[1]]]), t1$p.value))
cat(sprintf("%s: mean log2FC = %.3f, p-value = %.3e\n", lv[3], mean(log2FC[[lv[3]]]), t3$p.value))
cat(sprintf("%s: mean log2FC = %.3f, p-value = %.3e\n", lv[4], mean(log2FC[[lv[4]]]), t4$p.value))

cat("-----------------\n")
cat("RNA data\n")
cat("-----------------\n")

# now let's look at the RNA data
rna <- readRDS(rnaFile)

# calculate cell cycle scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# translate these to ensembl IDs
s.genes.ens <- mapIds(org.Hs.eg.db, keys=s.genes, column="ENSEMBL", keytype="SYMBOL")
g2m.genes.ens <- mapIds(org.Hs.eg.db, keys=g2m.genes, column="ENSEMBL", keytype="SYMBOL")
s.genes.ens <- s.genes.ens[!is.na(s.genes.ens)]
g2m.genes.ens <- g2m.genes.ens[!is.na(g2m.genes.ens)]

rna <- CellCycleScoring(rna, s.features = s.genes.ens, g2m.features = g2m.genes.ens, set.ident = TRUE)
# run SCTransform controlling for cell cycle effects
options(future.globals.maxSize= 100 * 1024^3)
t0 <- Sys.time()
rna <- SCTransform(rna, assay = "originalexp",
    vars.to.regress = c("S.Score", "G2M.Score"), 
    verbose = FALSE,
    conserve.memory = TRUE, return.only.var.genes = FALSE,
    ncells = 5000)
print(Sys.time()-t0)


browser()
saveRDS(rna, file=sprintf("%s/Sarropoulos_RNA_Seurat_cellCycleScored.rds", outDir))


DefaultAssay(rna) <- "SCT"
Idents(rna) <- rna$precisest_label
deg <- list()
for (i in seq_along(pairs)){
    cat(sprintf("* Finding marker genes for %s vs %s...\n", pairs[[i]][1], pairs[[i]][2]))
    deg[[i]] <- FindMarkers(
      object = rna,
      ident.1 = pairs[[i]][1],
      ident.2 = pairs[[i]][2],
      min.pct = 0.1,
      test.use = 'wilcox'
    )
}

dir_rna <- list()
log2FC_rna <- list()
for (i in seq_along(deg)){
    sig <- deg[[i]][deg[[i]]$p_val_adj < 0.05, ]
    nm <- paste(pairs[[i]], collapse="_vs_")
    log2FC_rna[[nm]] <- sig$avg_log2FC
    cat(sprintf("%s vs %s\n", pairs[[i]][1], pairs[[i]][2]))
    cat(sprintf("Num significant genes: %d\n", nrow(sig)))
    print(summary(sig$avg_log2FC))
    cat(sprintf("Num upregulated genes in group 1: %d\n", sum(sig$avg_log2FC > 0)))
    cat(sprintf("Num upregulated genes in group 2: %d\n", sum(sig$avg_log2FC < 0)))
    dir_rna[[i]] <- cbind(sum(sig$avg_log2FC > 0), sum(sig$avg_log2FC < 0))
    cat("--------\n")
}

# print a similar violin plot for the RNA log2 fold changes
log2FC_rna_df <- lapply(names(log2FC_rna), function(x) data.frame(pair=x, log2FC=log2FC_rna[[x]]))
log2FC_rna_df <- do.call(rbind, log2FC_rna_df)
log2FC_rna_df$pair <- factor(log2FC_rna_df$pair, levels=lv) 
p_rna <- ggplot(log2FC_rna_df, aes(x=pair, y=log2FC, fill=pair)) +
    geom_violin() + 
    geom_boxplot(width=0.1, outlier.shape=NA) +
    scale_fill_brewer(palette="Greens") +
    xlab("log2 fold change") +
    ggtitle("RNA: Distribution of log2 fold changes for significant genes") +
    theme_bw() +
    theme(axis.text=element_text(size=16), legend.position = "none")
ggsave(sprintf("%s/log2FC_rna_density.pdf", outDir), plot=p_rna, width=8, height=5)

# run WmW test for significantly lt 0 peaks in lv[1] and sig gt 0 for lv[3] and lv[4]
t1 <- t.test(log2FC_rna[[lv[1]]], mu = 0, alternative="greater")
t3 <- t.test(log2FC_rna[[lv[3]]], mu = 0, alternative="less")
t4 <- t.test(log2FC_rna[[lv[4]]], mu = 0, alternative="less")
cat("Significance tests for overall directionality of transcription\n")
cat(sprintf("%s: mean log2FC = %.3f, p-value = %.3e\n", lv[1], mean(log2FC_rna[[lv[1]]]), t1$p.value))
cat(sprintf("%s: mean log2FC = %.3f, p-value = %.3e\n", lv[3], mean(log2FC_rna[[lv[3]]]), t3$p.value))
cat(sprintf("%s: mean log2FC = %.3f, p-value = %.3e\n", lv[4], mean(log2FC_rna[[lv[4]]]), t4$p.value))

# plot a two-panel figure showing the directionality of changes in peaks and genes for each pair
# do this for the RNA layer as well, showing the number of upregulated genes in each group for each pair
dirChange_df <- as.data.frame(dirChange)
dirChange_df$pair <- rownames(dirChange_df)
dirChange_df_melt <- reshape2::melt(dirChange_df, id.vars="pair", variable.name="direction", value.name="num_peaks")
lv2 <- sub("_vs_", " vs ", lv)
dirChange_df_melt$pair <- factor(dirChange_df_melt$pair, levels=lv2)
# make scale log10 for the y-axis and add a small constant to avoid log(0)
# show each pair grouped
p_dirChange <- ggplot(dirChange_df_melt, aes(x=pair, y=num_peaks+1, fill=direction)) +
    geom_bar(stat="identity", position="dodge") +
    scale_y_log10() +
    scale_fill_manual(values=c("open in group 1"="steelblue", "open in group 2"="salmon")) +
    xlab("Comparison") +
    ylab("Number of significant peaks") +
    ggtitle("Directionality of significant peak changes") +
    theme_bw() +
    theme(axis.text=element_text(size=16), legend.position = "top") 
ggsave(sprintf("%s/directionality_peaks.pdf", outDir), plot=p_dirChange, width=12, height=5)

# now do the same for the RNA data
dirChange_rna_df <- as.data.frame(do.call(rbind, dir_rna))
colnames(dirChange_rna_df) <- c("upregulated in group 1", "upregulated in group 2")
dirChange_rna_df$pair <- lv2
dirChange_rna_df_melt <- reshape2::melt(dirChange_rna_df, id.vars="pair", variable.name="direction", value.name="num_genes")
dirChange_rna_df_melt$pair <- factor(dirChange_rna_df_melt$pair, levels=lv2)
p_dirChange_rna <- ggplot(dirChange_rna_df_melt, aes(x=pair, y=num_genes+1, fill=direction)) +
    geom_bar(stat="identity", position="dodge") +
    scale_y_log10() +
    scale_fill_manual(values=c("upregulated in group 1"="darkgreen", "upregulated in group 2"="lightgreen")) +
    xlab("Comparison") +
    ylab("Number of significant genes") +
    ggtitle("Directionality of significant gene changes") +
    theme_bw() +
    theme(axis.text=element_text(size=16), legend.position = "top")
ggsave(sprintf("%s/directionality_genes.pdf", outDir), plot=p_dirChange_rna, width=12, height=5)

# now print the number of cells in each group for each layer
cat("Num cells in each cluster:\n")
cat("ATAC\n")
cat("---------------\n")
cells_of_interest <- unique(unlist(pairs))
cat("Num cells in each cluster:\n")
tb <- table(atac$rna_precisest_label)
print(tb[cells_of_interest])
cat("\n\n")
cat("RNA\n")
cat("---------------\n")
tb <- table(rna$precisest_label)
print(tb[cells_of_interest])

}, error=function(ex){
    print(ex)
}, finally={
    sink()
})