rm(list=ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(glue)

default_wd <- this.path::this.dir()
setwd(default_wd) # set current scripts' dir as working dir
source("../../EM-seq/Analysis/FETHB2-FETHB3/utils.R")
source("../../EM-seq/Analysis/FETHB2-FETHB3/CpG/overlapEnrichment/utils.R")
source("./utils.R")

dt <- format(Sys.Date(),"%y%m%d")
dmr_date <- get_configs("CPG_DMR_DATE")
outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/{dmr_date}/oligoDesign")

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN

logFile <- sprintf("%s/additional_targets_%s.log",
                   outDir, 
                   dt
)


main <- function() {
  ### HARs overlapping fb CB H3K27ac & close to neuroDev/MB genes ###
  hars <- getHARs()
  names(hars) <- paste(seqnames(hars), start(hars), end(hars), sep = "_")
  fCB_h3k27ac <- getFetalCB_HistonePeaks()$H3K27ac_union
  
  
  hars_inter_fCBh3k27ac <- leftIntersect(hars, fCB_h3k27ac)
  hars_inter_fCBh3k27ac$name <- hars[hars_inter_fCBh3k27ac$origin]$name
  hars_inter_fCBh3k27ac <- getNearestGene(hars_inter_fCBh3k27ac, 
                                       gene_types = c("protein_coding")
  )
  
  # filter #
  g34_genes <- get_g34genes()
  neurodev_genes <- get_neurodevGenes()
  
  hars_inter_fCBh3k27ac$is_g34 <- hars_inter_fCBh3k27ac$nearestGene %in% g34_genes
  hars_inter_fCBh3k27ac$is_neurodev <- hars_inter_fCBh3k27ac$nearestGene %in% neurodev_genes
  
  res <- hars_inter_fCBh3k27ac[hars_inter_fCBh3k27ac$is_g34 | hars_inter_fCBh3k27ac$is_neurodev]
  names(res) <- paste(seqnames(res), start(res), end(res), sep = "_")
  
  large_regions <- res[width(res) > mpra_target_len]
  
  # trim #
  # get JASPAR TFBS of the large regions
  message("Getting JASPAR TFBS for the large regions")
  source("../../MISC/JASPAR_getTFBS/getTFBSMotifMatrix.R")
  jaspar_hg38 <- get_configs("JASPAR_DB")
  
  setwd("../../MISC/JASPAR_getTFBS/")
  tfbs_file <- sprintf("%s/additionalTargets_large_regions_jaspar_%s.out", outDir, dt)
  getTFBSMotifMatrix(large_regions, 
                     jaspar=jaspar_hg38, 
                     outFile = tfbs_file, 
                     tmpDir = outDir, 
                     convertMat = F
  )
  setwd(default_wd) 
  
  jaspar_tfbs <- read.table(tfbs_file)
  colnames(jaspar_tfbs) <- c("seqnames", "start", "end", 
                             "jaspar_id", "uncertain", "strand", "TF")
  jaspar_tfbs$strand <- NULL
  
  # keep only expressed in RL
  message("Removing TFBS of TFs likely not expressed in Hendrikse 2022 RL snRNA-seq")
  jaspar_tfbs <- makeGRangesFromDataFrame(jaspar_tfbs, keep.extra.columns = T)
  
  rl_genes <- rl_genes <- read.table(get_configs("RL_ACTIVE_GENES"))[[1]]
  jaspar_tfbs_rl <- jaspar_tfbs[jaspar_tfbs$TF %in% rl_genes]
  jaspar_tfbs_rl <- unique(jaspar_tfbs_rl)
  message(sprintf("%d unique RL TFBS in JASPAR", length(jaspar_tfbs_rl)))
  
  # for each large region, merge TFBS
  # only keep tfbs completely in the region
  region_tfbs_rl <- keep_ol(jaspar_tfbs_rl, 
                            large_regions, type = "within") 
  message(
    sprintf(
      "%d TFBS are completely within the large regions", 
      length(region_tfbs_rl)
    )
  )
  print(quantile(width(region_tfbs_rl)))
  
  olp <- findOverlapPairs(region_tfbs_rl, large_regions)
  tmp <- olp@first
  tmp$origin <- names(olp@second) # origin of larger_regions not dmr
  print(quantile(width(tmp)))
  
  splitted_tmp <- GenomicRanges::split(tmp, tmp$origin)
  
  # parallel stitch
  cl <- makeCluster(2)
  registerDoParallel(cl)
  
  ## cover all tfbs
  # select only one with first highest n_tfbs for each region
  reduced_tmp <- foreach(x = splitted_tmp, 
                         .combine = 'c', 
                         .export = c("stitchSelect_dmr_tfbs", "mpra_target_len")
  ) %dopar% {
    stitchSelect_dmr_tfbs(x, target_len = mpra_target_len)
  }
  
  stopCluster(cl)
  
  key <- reduced_tmp$origin
  reduced_tmp$name <- large_regions[key]$name
  reduced_tmp$GC <- large_regions[key]$GC
  reduced_tmp$nearestGene <- large_regions[key]$nearestGene
  reduced_tmp$is_g34 <- large_regions[key]$is_g34
  reduced_tmp$is_neurodev <- large_regions[key]$is_neurodev
  reduced_tmp$origin <- NULL

  # combine all regions #
  resized_smallRegions <- resize_interval(res[width(res) <= mpra_target_len], 
                              target_len = mpra_target_len)
  
  final <- c(resized_smallRegions, reduced_tmp)
  names(final) <- NULL
  
  
  outFile <- sprintf("%s/additionalTargets_%s.tsv", outDir, dt)
  write.table(as.data.frame(final), outFile, 
              col.names = T, row.names = T, quote = F, sep = "\t")
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
         finally = {
           message("\n\n--------- R sessionInfo ---------\n\n")
           print(sessionInfo())
           sink(type = "output")
           sink(type = "message")
           close(logFileCon)
         }
) 









