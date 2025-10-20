# call DMR using methylKit
library(methylKit)
rm(list=ls())

source("../../utils_PaiLab.R")

rootDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/ou
tput/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection"


dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/MethylKit_%s",rootDir, dt)
dir.create(outDir, recursive=FALSE)

cytoDir <- get_configs("CPG_REPORT_DIR")
#cytoDir <- "/home/rstudio/isilon/scratch/Rtmp1i7CTW"
phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv"

  ### Load data ###
  files <- list.files(cytoDir,, 
                      pattern = "snpFiltered.cytosine_report.txt$", 
                      full.names = T
  )

  sampleId <- extract_sampleId(files)
  tissueType <- get_TissueType(sampleId)
  
  cat("reading methylation data...\n")
  t0 <- Sys.time()
  methObj <- methRead(as.list(files),
                      sample.id= as.list(sampleId),
                      assembly = "hg38",
                      treatment = tissueType, 
                      dbtype="tabix",
                      context = "CpG",
                      pipeline = "bismarkCytosineReport", 
                      mincov = 5
  )
  
  print(Sys.time()-t0)

cat("Filtering by coverage...\n")
t0 <- Sys.time()
  filtered.methObj=filterByCoverage(methObj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)
print(Sys.time()-t0)

browser()
cat("Uniting methylation data...\n")
t0 <- Sys.time()
  methObj=unite(filtered.methObj, 
    destrand=TRUE, 
    min.per.group=4L,
  )
print(Sys.time()-t0)


browser()
###cat("Finding differentially methylated regions...\n")
###t0 <- Sys.time()
###  dmr <- calculateDiffMeth(methObj,
###                          mc.cores =4L)
###print(Sys.time()-t0)

cat("Extracting hyper/hypomethylated DMRs...\n")
t0 <- Sys.time()
hyper <- getMethylDiff(dmr, difference = 10, qvalue = 0.01, type = "hyper")
print(Sys.time()-t0)
cat("Number of hypermethylated DMRs:", nrow(hyper), "\n")
cat("now extracting hypomethylated DMRs...\n")
t0 <- Sys.time()
hypo <- getMethylDiff(dmr, difference = 10, qvalue = 0.01, type = "hypo")
print(Sys.time()-t0)
cat("Number of hypomethylated DMRs:", nrow(hypo), "\n")

