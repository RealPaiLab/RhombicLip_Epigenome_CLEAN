rm(list=ls())

library(methylKit)
library(ggplot2)
library(glue)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")

### dir config ###
rootDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/DMVs/CTsnv_excluded/%s/Segmentation",rootDir,dt)
logFile <- sprintf("%s/methylKit_segmentation_%s.log",outDir, dt)


main <- function() {
  message("MethylKit segmentation for each sample")
  message(sprintf("Output will be written under: %s\n",outDir))
  
  ### preprocess data ###
  filenames <- as.list(
    list.files(get_configs("CPG_REPORT_DIR"), 
               pattern="*.txt.gz$", 
               full.names=TRUE))
  sampleIds <- extract_sampleId(filenames)
  tissueType <- get_TissueType(unlist(sampleIds))
  
  message(sprintf("Loading data from \n- %s", paste(unlist(filenames), collapse = "\n- ")))
  methObj <- methRead(filenames,
                      sample.id= as.list(sampleIds),
                      assembly = "hg38",
                      treatment = tissueType, 
                      context = "CpG",
                      pipeline = "bismarkCytosineReport", 
                      mincov = 5
                      )
  
  ### plot segmentation diagnostic plots for one vz and one svz sample ###
  message("Testing segmentation diagnose for one SVZ sample.")
  for (i in methObj@.Data) {
    if (get_TissueType(i@sample.id) == 0) {
      test_svz <- i
      test_svz <- test_svz[test_svz$chr %in% paste0("chr", c(as.character(1:22), "X", "Y")),]
      break
    }
  }
  outFile <- sprintf("%s/testSVZ_segDiagnostic.png",outDir)
  png(outFile)
  res <- methSeg(test_svz, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
  dev.off()
  rm(test_svz)
  
  message("Testing segmentation diagnose for one VZ sample.")
  for (i in methObj@.Data) {
    if (get_TissueType(i@sample.id) == 1) {
      test_vz <- i
      test_vz <- test_vz[test_vz$chr %in% paste0("chr", c(as.character(1:22), "X", "Y")),]
      break
    }
  }
  outFile <- sprintf("%s/testVZ_segDiagnostic.png",outDir)
  png(outFile)
  res <- methSeg(test_vz, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
  dev.off()
  
  
  ### Formal segamentations ###
  message("Starting segmentation for each sample using 4 segmentation groups")
  for (i in methObj@.Data) {
    sample <- i@sample.id
    region <- ifelse(get_TissueType(sample) == 1, "VZ", "SVZ")
    
    t0 <- Sys.time()
    message(sprintf("Performing segmentation on sample %s, %s", sample, region))
    
    # segmentation #
    i <- i[i$chr %in% paste0("chr", c(as.character(1:22), "X", "Y")),]
    outFile <- sprintf("%s/%s_%s_segDiagnostic.png",outDir, sample, region)
    png(outFile)
    res <- methSeg(i,
                   diagnostic.plot=TRUE, 
                   maxInt=100, 
                   minSeg=10, 
                   G=1:4, 
                   join.neighbours = TRUE
                   ) # it seems minSeg is not working?
    dev.off()
    
    # plot seg group #
    outFile <- sprintf("%s/%s_%s_segScatter.png",outDir, sample, region)
    png(outFile)
    plot(res$seg.mean,
         log10(width(res)),pch=20,
         col=scales::alpha(rainbow(4)[as.numeric(res$seg.group)], 0.4),
         ylab="log10(width)",
         xlab="methylation proportion",
         main= sprintf("%s %s", sample, region)
    )
    dev.off()
    
    # convert to bed format #
    message("-- Converting segmentation to bed format")
    outFile <- sprintf("%s/%s_%s_seg.bed",outDir, sample, region)
    methSeg2bed(res,filename=outFile)
    
    
    message(Sys.time() - t0)
  }
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




