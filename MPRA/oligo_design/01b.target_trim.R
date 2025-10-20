rm(list=ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(glue)
library(foreach)
library(doParallel)

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

logFile <- sprintf("%s/target_trim_%s.log",
                   outDir, 
                   dt
)

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN



main <- function() {
  ### Data prep for target identification ###
  
  ## RL-VZ/SVZ DMRs ##
  dmr <- get_cpg_dmrs()
  names(dmr) <- paste(seqnames(dmr), start(dmr), end(dmr), sep = "_")
  dmr$source <- "dmr"

  ## dmr ol Aldinger human fetal hindbrain histone modifications (ChIP-seq) ##
  fetal_CB_enh <- get_fetalCB_enh(sprintf("%s/raw",
                                          get_configs("ALDINGER_CB_CRE_DIR"))
                                  )
  
  ## RL genes ##
  rl_genes <- read.table(get_configs("RL_ACTIVE_GENES"))[[1]]
  
  ## H3K27ac peak summit regions
  # to fit the target regions size for MPRA
  summit_gr <- get_fetalCB_h3k27ac_summit(
    up = floor(mpra_target_len/2), 
    down = mpra_target_len - floor(mpra_target_len/2)
    )
  
  
  ### Main ###
  ## Target region identification ##
  # intersect with Aldinger ChIP
  dmr_inter_fetalCbEnh <- leftIntersect(target_gr = dmr, ol_gr = fetal_CB_enh) 
  message(
    sprintf(
      "Found %d intersections (%i-%i, median %.2f) between DMRs and fetal CB enhancers in %d unique DMRs", 
      length(dmr_inter_fetalCbEnh), 
      min(width(dmr_inter_fetalCbEnh)),
      max(width(dmr_inter_fetalCbEnh)),
      median(width(dmr_inter_fetalCbEnh)),
      length(unique(dmr_inter_fetalCbEnh$origin))
    )
  )
  
  # get near summit region
  dmr_inter_fetalCbEnhSummit <- leftIntersect(target_gr = dmr_inter_fetalCbEnh, 
                                              ol_gr = summit_gr) 
  message(
    sprintf(
      "Found %d intersections (%i-%i, median %.2f) between DMRs and fetal CB enhancers and summits in %d unique DMRs", 
      length(dmr_inter_fetalCbEnhSummit), 
      min(width(dmr_inter_fetalCbEnhSummit)), 
      max(width(dmr_inter_fetalCbEnhSummit)),
      median(width(dmr_inter_fetalCbEnhSummit)),
      length(unique(dmr_inter_fetalCbEnhSummit$origin))
    )
  )
  
  # use Enhancer when there isn't summit overlap #
  fetalCbEnh_node <- drop_ol(dmr_inter_fetalCbEnh, dmr_inter_fetalCbEnhSummit)
  fetalCbEnh_node$source <- "dmr_fetalCbEnh"
  
  # use Summit when possible #
  summit_node <- dmr_inter_fetalCbEnhSummit
  summit_node$source <- "dmr_fetalCbEnhSummit"
  
  # further reduce regions larger than target size with TFBS #
  # replace dmrs trimmed 
  leftOver_dmr <- dmr[! names(dmr) %in% c(summit_node$origin, 
                                          fetalCbEnh_node$origin
                                          )
                      ]
  combined_regions <- c(leftOver_dmr, fetalCbEnh_node, summit_node)
  combined_regions[is.na(combined_regions$origin)]$origin <- names(
    combined_regions[is.na(combined_regions$origin)]
  )
  names(combined_regions) <- paste(seqnames(combined_regions),
                                   start(combined_regions),
                                   end(combined_regions), 
                                   sep = "_"
                                   )

  large_regions <- combined_regions[width(combined_regions) > mpra_target_len]
  names(large_regions) <- paste(seqnames(large_regions), 
                                start(large_regions), 
                                end(large_regions), 
                                sep = "_"
  )
  message(
    sprintf(
      "%d regions after combining all nodes; %d regions larger than %d bp", 
      length(combined_regions), 
      length(large_regions),
      mpra_target_len
    )
  )
  
  # get JASPAR TFBS of the large regions
  message("Getting JASPAR TFBS for the large regions")
  source("../../MISC/JASPAR_getTFBS/getTFBSMotifMatrix.R")
  jaspar_hg38 <- get_configs("JASPAR_DB")
  
  setwd("../../MISC/JASPAR_getTFBS/")
  tfbs_file <- sprintf("%s/large_regions_jaspar_%s.out", outDir, dt)
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
  cl <- makeCluster(16)
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
  reduced_tmp$origin <- large_regions[key]$origin
  reduced_tmp$source <- large_regions[key]$source
  reduced_tmp$source <- paste0(reduced_tmp$source, "_rlTFBS")
  
  # Final target regions #
  message("Generating final targets")
  targets <- c(combined_regions[width(combined_regions) <= mpra_target_len], 
               reduced_tmp
               )
  names(targets) <- NULL 
  targets$target <- paste(seqnames(targets), 
                          start(targets), 
                          end(targets), 
                          sep = "_"
  )
  print(quantile(width(targets)))
  message(sprintf("%i target regions after trimming", length(targets)))
  names(targets) <- NULL
  
  outFile <- sprintf("%s/trimmed_targets_%s.tsv", outDir, dt)
  write.table(as.data.frame(targets),
              outFile,
              col.names = T, 
              row.names = T,
              sep = "\t"
              )
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
