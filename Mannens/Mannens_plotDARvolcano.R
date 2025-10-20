# plot volcano plots of Mannens DARs

rm(list=ls())
library(Seurat)
library(ggplot2)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(AnnotationHub)

outRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024"
#inFile <- sprintf("%s/250612/Mannens_2024_seurat.qs", outRoot)

#darDir <- sprintf("%s/250619", outRoot)
#darDir <- sprintf("%s/250612", outRoot)
darDir <- sprintf("%s/250702", outRoot)

dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/DAR_volcano_plots", outRoot)
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

gpPairs <- list(
    g4=c("RL-VZ","RL-SVZ")
   # g2=c("RL-SVZ","Myeloid"),
    #g1=c("GCP","GC"),
    #g3=c("RL-SVZ","GCP")
    #g3=c("RL-SVZ","UBC precursors"),
)

for (cur in names(gpPairs)){
    g1 <- gpPairs[[cur]][1]
    g2 <- gpPairs[[cur]][2]
    print(sprintf("%s vs %s\n", g1, g2))

    curDir <- sprintf("%s/%s_vs_%s", outDir, g1, g2)
    da_peaks <- read.delim(sprintf("%s/dar_peaks_%s_vs_%s/da_peaks_%s_vs_%s.csv", 
        darDir,g1, g2, g1, g2),
        header=TRUE, sep="\t")
    browser()

      # plot volcano plot
  p <- ggplot(da_peaks, aes(x = avg_log2FC, y = -log10(p_val))) +
    geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.3), alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = sprintf("Volcano plot of DA peaks: %s vs %s", g1, g2),
         x = "Log2 Fold Change", y = "-log10 p-value")
  p <- p + theme(legend.position = "none")
pdf(sprintf("%s/volcano_%s_vs_%s.pdf", outDir, g1, g2), width = 6, height = 6)
print(p)
  dev.off()

  # plot the pvalue histogram 
  pdf(sprintf("%s/da_peaks_qqplot_%s_vs_%s.pdf", outDir, g1, g2), width = 6, height = 6)
  hist(da_peaks$p_val, breaks=200,
        main = sprintf("Histogram of p-values for DA peaks: %s vs %s", g1, g2),
        xlab = "p-value", ylab = "Frequency")
  dev.off()
}


#cat("Reading Mannens Seurat object from file\n")
#srat <- qs::qread(inFile)