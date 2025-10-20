rm(list=ls()); gc()

library(methylKit)
library(ape)
library(ggplot2)
library(glue)
library(ggrepel)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")

### Configs
dt <- format(Sys.Date(),"%y%m%d")

outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/EDA")
logFile <- glue("{outDir}/EDA_{dt}.log")


main <- function() {
  ### Load data ###
  files <- list.files(get_configs("CPG_REPORT_DIR"), 
                      pattern = "snpFiltered.cytosine_report.txt.gz$", 
                      full.names = T
  )
  sampleId <- extract_sampleId(files)
  tissueType <- get_TissueType(sampleId)
  
  methObj <- methRead(as.list(files),
                      sample.id= as.list(sampleId),
                      assembly = "hg38",
                      treatment = tissueType, 
                      context = "CpG",
                      pipeline = "bismarkCytosineReport", 
                      mincov = 1
  )
  
  
  ### Check quality ###
  for (i in 1:length(methObj)){
    print(methObj[[i]]@sample.id)
    getCoverageStats(methObj[[i]], both.strands=TRUE, plot=F)
    print("-----------------------------------------------------\n")
  }
  
  outFile <- glue("{outDir}/methylDistribution_{dt}.pdf")
  pdf(outFile)
  for (i in methObj) {
    getMethylationStats(i, plot=TRUE, both.strands=F)
  }
  dev.off()
  
  ### Merge ###
  # remove low count and extremely high count
  filtered.samples <- filterByCoverage(methObj, lo.count = 5, hi.perc = 99.5)
  rm(methObj); gc()
  filtered.norm.samples <- normalizeCoverage(filtered.samples, method = "median")
  meth <- unite(filtered.norm.samples, destrand=FALSE, mc.cores = 4)
  rm(filtered.norm.samples); gc()
  
  
  ### Clustering ###
  # get percent methylation matrix
  pm=percMethylation(meth)
  
  # calculate standard deviation of CpGs
  sds=matrixStats::rowSds(pm)
  
  # Visualize the distribution of the per-CpG standard deviation to determine a suitable cutoff
  outFile <- glue("{outDir}/sdsDistribution_{dt}.png")
  png(outFile)
  hist(sds, breaks = sqrt(length(sds)))
  dev.off()
  
  meth.filtered <- meth[sds > 50]
  hc <- clusterSamples(meth.filtered, dist="correlation", method="ward.D", plot = F)
  phylo <- as.phylo(hc)
  
  outFile <- glue("{outDir}/tree_methFiltered_{dt}.pdf")
  pdf(outFile)
  col <- c("orangered", "steelblue")[tissueType[phylo$tip.label] + 1]
  plot(phylo, label.offset = 0.02, tip.color = col, edge.width = 2, main = "meth.filtered")
  legend("topleft", legend = c("VZ", "SVZ"), fill = c("steelblue", "orangered"), 
         horiz = T, bg = "transparent", box.lwd = 0)
  dev.off()
  
  # PCA #
  pcs <- as.data.frame(PCASamples(meth.filtered, obj.return = T)$x)
  pcs$tissueType <- as.factor(tissueType[rownames(pcs)])
  pcs$sampleId <- rownames(pcs)
  
  .plot <- ggplot(pcs, aes(x = PC1, y = PC2, color = tissueType)) + 
    geom_point(alpha = 0.8) + 
    geom_text_repel(aes(label = sampleId), box.padding = 0.5)
  
  outFile <- glue("{outDir}/pca_{dt}.png")
  ggsave(outFile, .plot, dpi = 600, height = 6, width = 5)
  
  
  ### Check VZ/SVZ markers ###
  # SOX2 was used to delineate LCM boundry between VZ and SVZ 
  # TSS chr3:181711924-181711925 (hg38)
  # define promoter region -1000 ~ + 100
  
  chr <- "chr3"
  start_pos <- 181711925 - 10000
  end_pos <- 181711925 + 10000
  sox2_promoter <- meth[which(meth$chr == chr & 
                                meth$start >= start_pos & 
                                meth$end <= end_pos), 
  ]
  
  # Gene body
  # chr3:181711925-181714436
  # found standard deviation was zero for gene body
  print(nrow(sox2_promoter))
  
  sox2_hc <- clusterSamples(sox2_promoter, dist="correlation", method="ward.D", plot = F)
  sox2_phylo <- as.phylo(sox2_hc)
  
  outFile <- glue("{outDir}/tree_sox2_10kb_{dt}.pdf")
  pdf(outFile)
  col <- c("orangered", "steelblue")[tissueType[sox2_phylo$tip.label] + 1]
  plot(sox2_phylo, label.offset = 0.01, tip.color = col, edge.width = 2, type = "tidy", main = "SOX2 +/- 10kb")
  legend("topleft", legend = c("VZ", "SVZ"), fill = c("steelblue", "orangered"), 
         horiz = T, bg = "transparent", box.lwd = 0)
  dev.off()
  
  
  # EOMES TSS
  # chr3:27722322-27722323
  # chr3:27722497-27722498
  # chr3:27722710-27722711
  
  chr <- "chr3"
  start_pos <- 27722323 - 10000
  end_pos <- 27722711 + 10000
  eomes_promoter <- meth[which(meth$chr == chr & 
                                 meth$start >= start_pos & 
                                 meth$end <= end_pos), 
  ]
  print(nrow(eomes_promoter))
  
  
  eomes_hc <- clusterSamples(eomes_promoter, dist="correlation", method="ward.D", plot = F)
  eomes_phylo <- as.phylo(eomes_hc)
  
  outFile <- glue("{outDir}/tree_eomes_10kb_{dt}.pdf")
  pdf(outFile)
  col <- c("orangered", "steelblue")[tissueType[eomes_phylo$tip.label] + 1]
  plot(eomes_phylo, cex = 1.2, label.offset = 0.01, tip.color = col, edge.width = 2, type = "tidy", main = "EOMES +/- 10kb")
  legend("topleft", legend = c("VZ", "SVZ"), fill = c("steelblue", "orangered"), 
         horiz = T, bg = "transparent", box.lwd = 0)
  dev.off()
}




### main ###
if (! dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

logFileCon <- file(logFile, open = "wt")
sink(logFileCon, split = T, type = "output")
sink(logFileCon, type = "message")
tryCatch({main()}, 
         error = function(e) {message(e)},
         warning = function(w) {message(w)},
         finally = {
           message("\n\n--------- R sessionInfo ---------\n\n")
           print(sessionInfo())
           sink(type = "output")
           sink(type = "message")
           close(logFileCon)
         }
) 













