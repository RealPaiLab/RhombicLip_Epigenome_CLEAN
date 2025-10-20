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

logFile <- sprintf("%s/generate_targetProbes_%s.log",
                   outDir, 
                   dt
)

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN
budget <- configs$N_TARGET

additionalTargets_path <- sprintf("%s/additionalTargets_%s.tsv", 
                                  outDir, "240731")
rankedTargets_path <- sprintf("%s/targetRanking_%s.tsv", 
                              outDir, "240731")


### Utils ###
#' remove overlap with drop_gr from target_gr
remove_ol_region <- function(target_gr, drop_gr) {
  tmp <- lapply(split(target_gr), 
                function(x) {
                  # used df as do.call(c) couldn't work
                  as.data.frame(IRanges::setdiff(x, drop_gr))
                }
  )
  n <- names(tmp)[sapply(tmp, nrow) != 0]
  res <- GenomicRanges::makeGRangesFromDataFrame(do.call(rbind, tmp))
  names(res) <- n
  
  message(
    sprintf("Removed promoter overlapping targets: %d dropped because overlap",
            length(tmp) - length(res)
            )
    )
  
  return(res)
}


main <- function() {
  ### Load df ###
  ## additional target ##
  # process first because all need to be included
  additional_df <- read.table(additionalTargets_path, 
                              stringsAsFactors = F, header = T, sep = "\t"
  )
  rownames(additional_df) <- paste(additional_df$seqnames, 
                                   additional_df$start, additional_df$end,
                                   sep = "_"
                                   )
  
  
  ## ranked target ##
  targets_df <- read.table(rankedTargets_path, 
                           stringsAsFactors = F, header = T, sep = "\t"
  )
  rownames(targets_df) <- targets_df$target
  targets_df <- targets_df[! targets_df$ol_promoter,] # mpra probes not ol promoter
  
  
  ### Design tiled probes ###
  # not actually tiling because all target regions are <= mpra_target_len
  
  ## additional probes ##
  additional_targets <- GenomicRanges::makeGRangesFromDataFrame(additional_df)
  additionalProbe_list <- design_tiles(additional_targets, 
                                       width = mpra_target_len, overlap = 50)
  additionalProbes <- unlist(additionalProbe_list, use.names = F)
  
  # remove any OL with promoter
  additionalProbes <- drop_promoters(additionalProbes, 
                                     method = "overlap", 
                                     promoter_radius = 1e3)
  
  budget <- budget - length(additionalProbes)
  
  ## ranked probes ##
  targets <- GenomicRanges::makeGRangesFromDataFrame(
    targets_df[1:min(budget+20, nrow(targets_df)),],
    keep.extra.columns = T
  )

  # will deal with USH2A manually later
  targets <- targets[! (targets$nearestGene == "USH2A" & targets$ol_har)]
  
  targetProbe_list <- design_tiles(targets, 
                                   width = mpra_target_len, overlap = 50)
  targetProbes <- unlist(targetProbe_list, use.names = F)
  
  targetProbes <- drop_promoters(targetProbes, 
                                 method = "overlap", 
                                 promoter_radius = 1e3)
  
  targetProbes <- targetProbes[1:(budget-2)] # leave 2 budgets for USH2A wild/mut
  
  ## summarize ##
  message(sprintf("There are %d overlaps between targetProbes and additionalProbes",
                  length(findOverlaps(targetProbes, additionalProbes))
  )
  )
  
  message(
    sprintf(
      "%d overlaps within targetProbes", 
      length(targetProbes) - length(reduce(targetProbes))
    )
  )
  
  message(
    sprintf(
      "%d overlaps within additionalProbes",
      length(additionalProbes) - length(reduce(additionalProbes))
    )
  )
  
  ### Combine probes ###
  colnames(additional_df)[colnames(additional_df) == "is_g34"] <- "nearestGeneIsMbGene"
  colnames(additional_df)[colnames(additional_df) == "is_neurodev"] <- "nearestGeneIsNeurodevGene"
  additional_meta_cols <- c("nearestGene", "name", 
                            "nearestGeneIsMbGene", "nearestGeneIsNeurodevGene", 
                            "n_TFBS")
  additionalProbes$class <- "additional"
  for (c in additional_meta_cols) {
    elementMetadata(additionalProbes)[[c]] <- additional_df[names(additionalProbes), c]
  }
  
  target_meta_cols <- colnames(targets_df)[c(8, 9, 10, 11, 15:ncol(targets_df))]
  targetProbes$class <- "ranked"
  for (c in target_meta_cols) {
    elementMetadata(targetProbes)[[c]] <- targets_df[names(targetProbes), c]
  }
  
  
  # add USH2A wild type
  ush2a_mutPos <- GenomicRanges::GRanges(
    seqnames = "chr1", 
    ranges = IRanges(start = 216183532, end = 216183533)
    )
  ush2a_targetRegion <- resize_interval(ush2a_mutPos, 
                                        target_len = mpra_target_len
                                        )
  names(ush2a_targetRegion) <- "USH2A_WT_G"
  ush2a_targetRegion$class <- "SNV"
  
  combinedProbes <- c(additionalProbes, targetProbes, ush2a_targetRegion)
  
  message(
    sprintf(
      "Selected %d probes targeting %d intervals: Addtional - %d,%d; Ranked - %d,%d",
      length(combinedProbes),
      length(unique(names(combinedProbes))),
      table(combinedProbes$class)["additional"],
      length(unique(names(combinedProbes[combinedProbes$class == "additional"]))),
      table(combinedProbes$class)["ranked"],
      length(unique(names(combinedProbes[combinedProbes$class == "ranked"])))
    )
  )
  message(
    sprintf("Probe length: Min (%d), Median (%d), max (%d)",
            min(width(combinedProbes)),
            median(width(combinedProbes)),
            max(width(combinedProbes))
    )
  )
  
  
  ## output GRanges df ##
  outFile <- sprintf("%s/targetProbes_hg38_%s.rds", outDir, dt)
  saveRDS(combinedProbes, outFile)
  
  ## output bed ##
  outFile <- sprintf("%s/targetProbes_hg38_%s.bed", outDir, dt)
  elementMetadata(combinedProbes) <- NULL
  export.bed(combinedProbes, outFile)
  
  ## output fasta ##
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  seqs <- getSeq(genome, combinedProbes)
  
  # !ADD USH2A wild/mut seq #
  # target mut is chr1:216183532-216183533 G>A
  # design target region center on the mut pos
  ush2a_wt_seq <- getSeq(genome, ush2a_targetRegion)
  relative_loc <- 216183533 - start(ush2a_targetRegion) + 1
  ush2a_mut_seq <- Biostrings::replaceAt(ush2a_wt_seq, 
                             at = IRanges(relative_loc, relative_loc), 
                             "A"
                             )
  names(ush2a_mut_seq) <- "USH2A_MUT_G2A"
  
  # combine all
  seqs <- c(seqs, ush2a_mut_seq)
  
  # check restriction sites
  check_SceI(seqs)
  
  outFile <- sprintf("%s/targetProbes_hg38_%s.fasta", outDir, dt)
  memes::write_fasta(seqs, outFile)
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

