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

logFile <- sprintf("%s/target_selection_%s.log",
                   outDir, 
                   dt
)

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN


main <- function() {
  ### Data prep for target identification ###
  
  ## RL-VZ/SVZ DMRs ##
  dmr <- get_cpg_dmrs()

  ## dmr ol Aldinger human fetal hindbrain histone modifications (ChIP-seq) ##
  fetal_CB_enh <- get_fetalCB_enh()
  
  ## Smith et al G3/4 MB ATAC-seq ##
  g3_atac <- import.bed("/data/xsun/20240314/cre/Smith2022-Group3/allEnhancerLikeElements.bed")
  g4_atac <- import.bed("/data/xsun/20240314/cre/Smith2022-Group4/allEnhancerLikeElements.bed")
  
  ## RL genes ##
  rl_genes <- read.table("/data/xsun/db/meme/HOCOMOCOv12/activeGenes/Hendrikse2022_RL_activeGenes")[[1]]
  
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
  fetalCbEnh_node$source <- "fetalCbEnh"
  
  # use Summit when possible #
  summit_node <- dmr_inter_fetalCbEnhSummit
  summit_node$source <- "summit"
  
  # other non-overlapping DMRs use tumour atac-seq #
  # get DMRs having fetal enhancer
  matched_dmr <- GenomicRanges::makeGRangesFromDataFrame(
    do.call(rbind, 
            lapply(strsplit(unique(dmr_inter_fetalCbEnh$origin), "_"), 
                   function(x) {as.data.frame(t(x))})
    ) %>% rename(seqnames = V1, start = V2, end = V3)
  )
  # exclude DMRs
  orphan_dmr <- GenomicRanges::setdiff(dmr, matched_dmr)
  message(
    sprintf(
      "%d out of %d DMRs do not overlap with fetal hindbrain enhancers.", 
      length(orphan_dmr), 
      length(dmr)
    )
  )
  
  # intersect with G3/4 ATAC
  atac_node <- leftIntersect(orphan_dmr, 
                             GenomicRanges::union(g3_atac, g4_atac))
  
  message(
    sprintf(
      "Found %d intersections (%i-%i, median %.2f) between orphan DMRs and G3/4 MB ATAC peaks in %d unique DMRs", 
      length(atac_node), 
      min(width(atac_node)),
      max(width(atac_node)),
      median(width(atac_node)),
      length(unique(atac_node$origin))
    )
  )
  atac_node$source <- "atac"
  
  # Check dmrs still have no match #
  tmp <- c(fetalCbEnh_node, summit_node, atac_node)
  tmp_matched_dmr <- GenomicRanges::makeGRangesFromDataFrame(
    do.call(rbind, 
            lapply(strsplit(unique(tmp$origin), "_"),
                   function(x) {as.data.frame(t(x))})
    ) %>% rename(seqnames = V1, start = V2, end = V3)
  )
  leftOver_dmr <- GenomicRanges::setdiff(dmr, tmp_matched_dmr)
  leftOver_dmr$origin <- paste(seqnames(leftOver_dmr), 
                               start(leftOver_dmr), 
                               end(leftOver_dmr), 
                               sep = "_")
  leftOver_dmr$source <- "leftOver"
  
  # further reduce regions larger than target size with TFBS # 
  combined_regions <- c(fetalCbEnh_node, summit_node, atac_node, leftOver_dmr)
  
  large_regions <- combined_regions[width(combined_regions) > mpra_target_len]
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
  jaspar_hg38 <- "/.mounts/labs/pailab/src/ucsc-goldenpath/hg38/JASPAR2024.bb"
  
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
  max_tfbs_size <- max(width(region_tfbs_rl))
  message(sprintf("Max size of tfbs is %d", max_tfbs_size))
  
  large_regions_cp <- large_regions
  # remove origin for leftIntersect to track; 
  # CANDO: change function to ignore when needed
  large_regions_cp$origin <- NULL 
  tmp <- leftIntersect(large_regions_cp, region_tfbs_rl)
  print(quantile(width(tmp)))
  
  splitted_tmp <- GenomicRanges::split(tmp, tmp$origin)
  reduced_splitted_tmp <- lapply(
    splitted_tmp, 
    function(x) {
      GenomicRanges::reduce(x, 
                            min.gapwidth= mpra_target_len - max_tfbs_size*2 + 1
                            )
    }
  )
  ln <- unlist(lapply(reduced_splitted_tmp, length))
  print(table(ln))
  
  reduced_tmp <- unlist(as(reduced_splitted_tmp, "GRangesList"))
  names(large_regions) <- paste(seqnames(large_regions), 
                                start(large_regions), 
                                end(large_regions), 
                                sep = "_"
  )
  reduced_tmp$origin <- large_regions[names(reduced_tmp)]$origin
  reduced_tmp$source <- paste0(large_regions[names(reduced_tmp)]$source, 
                               "_rlTFBS"
  )
  
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
  message(sprintf("%i target regions were selected.", nrow(targets)))
  
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
