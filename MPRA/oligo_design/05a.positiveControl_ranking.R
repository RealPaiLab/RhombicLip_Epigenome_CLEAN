rm(list=ls())

library(glue)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(msigdbr)
library(ggplot2)
library(ggpubr)

default_wd <- this.path::this.dir()
setwd(default_wd) # set current scripts' dir as working dir
source("../../EM-seq/Analysis/FETHB2-FETHB3/CpG/overlapEnrichment/utils.R")
source("../../EM-seq/Analysis/FETHB2-FETHB3/utils.R")
source("./utils.R")

dt <- format(Sys.Date(),"%y%m%d")
dmr_date <- get_configs("CPG_DMR_DATE")

outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/{dmr_date}/oligoDesign")

logFile <- sprintf("%s/posCtrl_ranking_%s.log",
                   outDir,
                   dt
)

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN
budget <- configs$N_POS

crs_bed <- sprintf("%s/targetProbes_hg38_%s.bed", outDir, 240731)

### Functions ###
loadNarrowPeaks <- function(narrowPeakFiles, builds = NA) {
  message(sprintf("%i narrowPeak files provided.", length(narrowPeakFiles)))
  loaded <- mapply(
    function(x, y) {
      message(glue("--- Importing from {x}."))
      res <- rtracklayer::import(x)
      
      if (! is.na(y)) {
        if (y == "hg19") {
          res <- liftOver_gr(res)
        }
      }
      
      return(res)
      },
    narrowPeakFiles, 
    builds
    )
  
  return(loaded)
  }


#' Find consensus peaks across multiple narrowPeak files 
#' @param narrowPeakFiles (list) A list of paths to narrowPeak files
#' @param minCov (numeric) The minimum coverage to be counted as consensus. Default
#' covering all files
#' @return (GRanges) The consensus peaks
findConsensusPeaks <- function(narrowPeakList, minCov = NA) {
  # preprocess #
  tmp <- lapply(narrowPeakList, 
                function(x) {
                  len_original_x <- length(x)
                  # reduce to ensure homogeneous max coverage
                  x <- IRanges::reduce(x)
                  message(sprintf("%d intervals (%d before reduction) with a median width of %.2f",
                                  length(x),
                                  len_original_x,
                                  median(width(x))
                                  )
                  )
                  
                  return(x)
                  }
                )
  
  
  # set coverage #
  if (is.na(minCov)) {
    minCov <- length(narrowPeakList)
  }
  
  if (minCov > length(narrowPeakList)) {
    stop("minCov cannot be larger than the number of narrowPeak files.")
  }
  
  message(glue("Finding ranges with >= {minCov} coverage."))
  
  # find consensus #
  narrowPeakList <- GRangesList(tmp)
  peakCov <- IRanges::coverage(narrowPeakList)
  covered_ranges <- IRanges::slice(peakCov, lower=minCov, rangesOnly=T)
  res <- GRanges(covered_ranges)
  message(sprintf("Found %i consensus peaks with %.2f median width", 
                  length(res),
                  median(width(res))
                  )
          )

  return(res)
}


#' Generate a matrix of consensus peaks SignalValue (not real as peak may be different from original)
#' @param narrowPeakList 
#' @param consensusPeaks (GRanges) The consensus peaks
#' @param (matrix) Row: consensus peak intervals. Column: narrowPeak file/sample
get_PeakSignalValue_matrix <- function(narrowPeakList, consensusPeaks) {
  # Create base mat #
  mat <- matrix(nrow = length(consensusPeaks), 
                ncol = length(narrowPeakList)
  )
  rownames(mat) <- paste(seqnames(consensusPeaks), 
                         start(consensusPeaks), 
                         end(consensusPeaks), 
                         sep = "_"
  )
  colnames(mat) <- names(narrowPeakList)
  
  # loop over all files #
  for (filename in names(narrowPeakList)) {
    narrowPeak <- narrowPeakList[[filename]]
    # e.g. ipsc_atac contains a/b of the same peak, only keep the first one
    narrowPeak <- IRanges::unique(narrowPeak)
    # min-max scaling as signalValues may be calculated differently by authors
    m <- min(narrowPeak$signalValue)
    d <- max(narrowPeak$signalValue) - m
    narrowPeak$signalValueScaled <- (narrowPeak$signalValue - m)/d
    
    ol <- findOverlaps(consensusPeaks, narrowPeak)
    # the order won't change as consusPeaks is the query and all are overlapped
    mat[, filename] <- narrowPeak[subjectHits(ol)]$signalValueScaled
  }
  
  return(mat)
}


#' Annotate intervals and rank the intervals based on PeakSignalValue 
#' @param mat (matrix) The PeakSignalValue matrix, rownames peak, colnames sample
#' @param max_size (numeric) maximum size of intervals will be kept. Default 171 bp for MPRA
#' @return (data.frame)
rank_intervals <- function(mat, max_size = 171) {
  # rank regions
  mat_scaled <- apply(mat, 2, function(x) {
    (x-min(x, na.rm = T))/((max(x, na.rm = T) - min(x, na.rm = T)))
  })
  scaledSignalValue_mean <- rowMeans(mat_scaled)

  splitted_names <- stringr::str_split(rownames(mat_scaled), 
                                       pattern = "_", 
                                       simplify = T
  )
  df <- data.frame(seqnames = splitted_names[,1],
                   start = as.integer(splitted_names[,2]),
                   end = as.integer(splitted_names[,3]))
  rownames(df) <- rownames(mat_scaled)
  message("--- Getting nearest genes for each interval")
  df_gr <- GenomicRanges::makeGRangesFromDataFrame(df)
  
  df_gr <- getNearestGene(df_gr, gene_types = c("protein_coding"))
  df_gr$scaledSignalValue_mean <- scaledSignalValue_mean[names(df_gr)]
  df_gr$isHousekeeping <- df_gr$nearestGene %in% get_housekeepingGenes()
  df_gr$isNeurodev <- df_gr$nearestGene %in% get_neurodevGenes()
  df_gr$isNH <- df_gr$isHousekeeping | df_gr$isNeurodev
  
  res <- as.data.frame(df_gr) %>%
    mutate(width = width(df_gr)) %>%
    filter(width <= max_size) %>%
    arrange(
      desc(scaledSignalValue_mean),
      desc(isHousekeeping),
      desc(isNeurodev)
    ) %>%
    mutate(rank = row_number())
  
  message(sprintf("%i intervals <= %i bp", nrow(res), max_size))
  message(sprintf("%i are near housekeeping genes, %i neurodevelopment genes",
                  sum(res$isHousekeeping),
                  sum(res$isNeurodev)
  ) 
  )
  
  return(res)
}


## Functions for Whalen dataset processing
#' 
get_whalen_paths <- function(dir) {
  files <- list.files(
    dir, 
    "^GSM30160[0-9]{2}_human-(WTc|HS1-11)-N[2,3]-[1-3]-(DNA|RNA).tsv.gz"
    )
  
  exps <- sort(unique(stringr::str_extract(files, "human-(WTc|HS1-11)-N[2,3]-[1-3]")))
  message(sprintf("Found %d pair of DNA/RNA results from Whalen 2023 MPRA:", 
                  length(exps)
                  ))
  message(paste(exps, collapse = "\n"))
  
  res <- list()
  for (exp in exps) {
    res[[exp]][["dna"]] <- sprintf("%s/%s",
                                   dir,
                                   files[grepl(glue("{exp}-DNA"), files)]
                                   )
    res[[exp]][["rna"]] <- sprintf("%s/%s",
                                   dir,
                                   files[grepl(glue("{exp}-RNA"), files)]
                                   )
  }
  
  return(res)
}

#' load tsv files from Whalen 2023 paper and covert to wide matrix
whalen_long2wide <- function(path) {
  d <- read.table(path, stringsAsFactors = F, header = F)
  colnames(d) <- c("obs", "count", "seqname")
  
  tmp <- stringr::str_split(d$seqname, "__:", simplify = T)
  d$seq <- tmp[,1]
  d$num <- tmp[,2]
  
  d_wide <- reshape2::dcast(d, seq ~ num, value.var = "count", fill = 0)
  rownames(d_wide) <- d_wide$seq
  d_wide$seq <- NULL
  
  return(as.matrix(d_wide))
}


#' normalize and get rna/dna ratio 
#' @param path_list (character) A list structured as {sample:{dna:x, rna:y}}
#' @return (character) A list of samples' normalized dna, rna, ratio, and log2
prepare_whalen <- function(path_list) {
  res <- list()
  for (i in names(path_list)) {
    dna <- whalen_long2wide(path_list[[i]][["dna"]])
    # adopted from MPRAflow merge_label.py
    dna <- (rowSums(dna) + 1)/(rowSums(dna != 0) + 1)/sum(dna)*1e6
    rna <- whalen_long2wide(path_list[[i]][["rna"]])
    rna <- (rowSums(rna) + 1)/(rowSums(rna != 0) + 1)/sum(rna)*1e6
    
    res[[i]][["dna"]] <- dna
    res[[i]][["rna"]] <- rna
    res[[i]][["ratio"]] <- rna/dna
    res[[i]][["log2"]] <- log2(rna/dna)
  }
  
  return(res)
}

#'
merge_whalen_ratios <- function(input) {
  tmp <- data.frame(matrix(NA,
                           nrow = max(unlist(lapply(input, 
                                                    function(x) {length(x$dna)}
                           )
                           )
                           ), 
                           ncol = length(input)
  ),
  row.names = unique(unlist(lapply(input,
                                   function(x) {names(x$dna)}
  )))
  )
  colnames(tmp) <- names(input)
  for (sample in names(input)) {
    tmp[,sample] <- input[[sample]][["ratio"]][rownames(tmp)]
  }
  
  tmp$class <- tidyr::replace_na(stringr::str_extract(rownames(tmp), 
                                                      "(POSITIVE|NEGATIVE)"),
                                 "CRS"
                                 )
  tmp$isENCODE <- grepl("ENCODE", rownames(tmp))
  
  tmp$mean <- rowMeans(tmp[,1:length(input)], na.rm = T)
  
  return(tmp)
}

#' plot positive and negative controls rna/dna ratio
plot_whalen_posNeg <- function(df, quantile = c(0.1, 0.9)) {
  .plot <- ggplot(df[df$class != "CRS",], aes(y = mean, x = class)) + 
    geom_violin(draw_quantiles = quantile) + 
    geom_jitter(aes(color = isENCODE), 
                position = position_jitter(seed = 1234, width = 0.1)) +
    scale_color_manual(values = c("TRUE" = alpha("red", 0.6), 
                                  "FALSE" = alpha("black", 0.2))) + 
    stat_compare_means(method = "wilcox.test") + 
    theme_bw()
  
  return(.plot)
}



main <- function() {
  ### set up peak input paths ###
  ## H3K27ac
  # hg19
  N2_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                          "/ChIP-seq/Whalen_2023",
                          "/GSE110758_human-HS1-11-N2-pooled-K27ac.narrowPeak.gz")
  # hg19
  N3_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                          "/ChIP-seq/Whalen_2023",
                          "/GSE110758_human-HS1-11-N3-pooled-K27ac.narrowPeak.gz")
  
  # hg38
  wct11_ipsc_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                                  "/ChIP-seq/Ren_2021",
                                  "/GSE166835_H3K27ac_peaks.narrowPeak.gz")
  
  # hg19
  npc_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                           "/ChIP-seq/Choi_2021",
                           "/GSE158382_NPC_h3k27ac_merged_sorted_noDup_peaks.narrowPeak.gz")
  # hg19
  ipsc_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                            "/ChIP-seq/Choi_2021",
                            "/GSE158382_iPSC_h3k27ac_merged_sorted_noDup_peaks.narrowPeak.gz")
  
  ## ATAC ##
  # hg19
  ipsc_atac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                         "/ATAC-seq/Choi_2021",
                         "/GSE158382_iPSC_ATAC_peaks.narrowPeak.gz")
  # hg19
  npc_atac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                        "/ATAC-seq/Choi_2021",
                        "/GSE158382_NPC_ATAC_peaks.narrowPeak.gz")
  
  ## path list ##
  narrowPeakFiles <- list(N2_h3k27ac = N2_h3k27ac_file,
                          N3_h3k27ac = N3_h3k27ac_file,
                          npc_h3k27ac = npc_h3k27ac_file,
                          npc_atac = npc_atac_file
  )
  
  
  ### Consensus intervals ###
  narrowPeakList <- loadNarrowPeaks(narrowPeakFiles, 
                                    builds = c("hg19", "hg19", "hg19", "hg19")
  )
  
  narrowPeakList_tmp <- narrowPeakList
  narrowPeakList_tmp$feCB_h3k27ac <- getFetalCB_HistonePeaks()$H3K27ac
  encode_cres <- get_encode_cres()
  narrowPeakList_tmp$ENCODE_ELS <- encode_cres[grepl("ELS", encode_cres$Type)]
  
  consensusPeaks <- findConsensusPeaks(narrowPeakList_tmp)
  consensusPeaks <- drop_promoters(consensusPeaks, 
                                   method = "overlap", 
                                   promoter_radius = 1000
                                   )
  
  mat <- get_PeakSignalValue_matrix(narrowPeakList, # no signalValue in els/feCB
                                    consensusPeaks
  )
  intervalRank_df <- rank_intervals(mat, max_size = mpra_target_len)
  # remove those too small
  intervalRank_df <- intervalRank_df[intervalRank_df$width > 10,]
  
  ### VISTA hindbrain enhancer ###
  all_lines <- readLines(get_configs("VISTA_HB"))
  filtered_lines <- grep("^>", all_lines, value = TRUE)
  loc <- stringr::str_split(filtered_lines, "\\|", simplify = T)[,2]
  loc <- stringr::str_trim(loc)
  element <- stringr::str_split(filtered_lines, "\\|", simplify = T)[,3]
  element <- sub(" ", "_", stringr::str_trim(element))
  element <- glue("VISTA_{element}")
  
  
  vista_hb <- as.data.frame(stringr::str_split(loc, ":|-", simplify = T))
  rownames(vista_hb) <- element
  colnames(vista_hb) <- c("seqnames", "start", "end")
  vista_hb_gr_hg19 <- GenomicRanges::makeGRangesFromDataFrame(vista_hb)
  vista_hb_gr_hg38 <- liftOver_gr(vista_hb_gr_hg19)
  
  ## trim vista ##
  # all vista hindbrain enhancers larger than 500bp
  # get JASPAR TFBS of the large regions
  message("Getting JASPAR TFBS for the large regions")
  source("../../MISC/JASPAR_getTFBS/getTFBSMotifMatrix.R")
  jaspar_hg38 <- get_configs("JASPAR_DB")
  
  setwd("../../MISC/JASPAR_getTFBS/")
  tfbs_file <- sprintf("%s/vista_large_regions_jaspar_%s.out", outDir, dt)
  getTFBSMotifMatrix(vista_hb_gr_hg38, 
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
  
  rl_genes <- read.table(get_configs("RL_ACTIVE_GENES"))[[1]]
  jaspar_tfbs_rl <- jaspar_tfbs[jaspar_tfbs$TF %in% rl_genes]
  jaspar_tfbs_rl <- unique(jaspar_tfbs_rl)
  message(sprintf("%d unique RL TFBS in JASPAR", length(jaspar_tfbs_rl)))
  
  # for each large region, merge TFBS
  # only keep tfbs completely in the region
  region_tfbs_rl <- keep_ol(jaspar_tfbs_rl, 
                            vista_hb_gr_hg38, type = "within") 
  message(
    sprintf(
      "%d TFBS are completely within the large regions", 
      length(region_tfbs_rl)
    )
  )
  print(quantile(width(region_tfbs_rl)))
  
  olp <- findOverlapPairs(region_tfbs_rl, vista_hb_gr_hg38)
  tmp <- olp@first
  tmp$origin <- names(olp@second) # origin of larger_regions not dmr
  print(quantile(width(tmp)))
  
  splitted_tmp <- GenomicRanges::split(tmp, tmp$origin)
  
  # parallel stitch
  cl <- makeCluster(10)
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
  
  tmp <- reduced_tmp
  names(tmp) <- tmp$origin
  elementMetadata(tmp) <- NULL
  vista_df <- as.data.frame(tmp)
  vista_df[,colnames(intervalRank_df)[6:ncol(intervalRank_df)]] <- NA
  vista_df$rank <- 0
  

  ### get Whalen et al 2023 controls
  whalen_files <- get_whalen_paths(get_configs("WHALEN_MPRA_DIR"))
  whalen_processed <- prepare_whalen(whalen_files)
  whalen_ratio_df <- merge_whalen_ratios(whalen_processed)
  
  # select consistent controls #
  # neg
  neg <- whalen_ratio_df[whalen_ratio_df$class == "NEGATIVE",]
  negCtrl <- rownames(neg[neg$mean < quantile(neg$mean, 0.1) & neg$isENCODE,])
  
  negCtrl_tmp <- stringr::str_split(negCtrl, "___", simplify = T)[,1]
  negCtrl_bed <- as.data.frame(stringr::str_split(negCtrl_tmp, 
                                                  "[:|-]", 
                                                  simplify = T)
                               )
  colnames(negCtrl_bed) <- c("seqnames", "start", "end")
  negCtrl_bed$name <- negCtrl
  
  outFile <- sprintf("%s/whalen_MPRA_selected_negCtrl_hg19_%s.bed", outDir, dt)
  write.table(negCtrl_bed, outFile, 
              row.names = F, col.names = F, 
              quote = F, sep = "\t")
  
  # whalen pos
  pos <- whalen_ratio_df[whalen_ratio_df$class == "POSITIVE",]
  posCtrl <- rownames(pos[pos$mean > quantile(pos$mean, 0.9) & pos$isENCODE,])
  
  posCtrl_tmp <- stringr::str_split(posCtrl, "___", simplify = T)[,1]
  posCtrl_bed <- as.data.frame(stringr::str_split(posCtrl_tmp, 
                                                  "[:|-]", 
                                                  simplify = T)
                               )
  colnames(posCtrl_bed) <- c("seqnames", "start", "end")
  posCtrl_bed$name <- posCtrl
  
  outFile <- sprintf("%s/whalen_MPRA_selected_posCtrl_hg19_%s.bed", outDir, dt)
  write.table(posCtrl_bed, outFile, 
              row.names = F, col.names = F, 
              quote = F, sep = "\t")
  
  posCtrl_gr <- import.bed(outFile) # for downstream combination
  posCtrl_gr <- liftOver_gr(posCtrl_gr)
  
  # save plot 
  .plot <- plot_whalen_posNeg(whalen_ratio_df)
  outFile <- sprintf("%s/whalen_MPRA_control_ratio_%s.png", outDir, dt)
  ggsave(outFile, .plot, width = 7, height = 6, dpi = 600)

  
  ### Combine posCtrl and ranked intervals
  posCtrl_df <- as.data.frame(posCtrl_gr)
  posCtrl_df[,colnames(intervalRank_df)[6:ncol(intervalRank_df)]] <- NA
  posCtrl_df$rank <- 0
  rownames(posCtrl_df) <- posCtrl_df$name
  posCtrl_df$name <- NULL
  
  ### Combine all ###
  final_df <- rbind(posCtrl_df, intervalRank_df, vista_df)
  final_df <- arrange(final_df, rank)
  
  ### Resize intervals ###
  resized_df <- resize_interval(final_df[,1:3], mpra_target_len)
  colnames(resized_df) <- paste0("resized_", colnames(resized_df))
  final <- cbind(final_df,
                 resized_df
  )
  
  # write rankings #
  outFile <- sprintf("%s/posCtrlRanking_%s.tsv", outDir, dt)
  write.table(final, outFile, col.names = T, row.names = F, quote = F, sep = "\t")
  
  # write fasta seq #
  message("removing problematic intervals")
  resized_final <- final %>%
    dplyr::select(resized_seqnames, resized_start, resized_end) %>%
    rename(seqnames = resized_seqnames, 
           start = resized_start,
           end = resized_end
           )
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(resized_final)
  gr <- unique(gr) # remove duplicated intervals, keep only the first
  
  # remove positive controls that are exactly the same as CRS
  crs <- import.bed(crs_bed)
  exact_matches <- names(findOverlapPairs(gr, crs, type = "equal")@first)
  gr <- gr[! names(gr) %in% exact_matches]
  
  # extract seqs
  message("Extracting sequences from hg38 genome")
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  seqs <- getSeq(genome, gr[1:budget])
  
  check_SceI(seqs)
  
  outFile <- sprintf("%s/posCtrlRanking_hg38_%s.fasta", outDir, dt)
  message(sprintf("Writing sequences to %s", outFile))
  memes::write_fasta(seqs, outFile)
  
  # write bed #
  outFile <- sprintf("%s/posCtrlRanking_hg38_%s.bed", outDir, dt)
  export.bed(gr[1:budget], outFile)
}



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




