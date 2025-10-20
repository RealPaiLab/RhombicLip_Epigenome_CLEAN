rm(list=ls())

library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rtracklayer)
library(glue)
library(dplyr)

setwd(this.path::this.dir())
source("../utils.R")
source("../../../EM-seq/Analysis/FETHB2-FETHB3/utils.R")

dt <- format(Sys.Date(),"%y%m%d")

### Configs ###
outDir <- "/.mounts/downloads/private/xinghansun/UCSC_track_hub/humanFetalRhombicLip"
tFile <- sprintf("%s/trackDb.txt", outDir)
#system(sprintf("cat /dev/null > %s",tFile)) # !! don't run unless certain
srcUrl <- "https://downloads.res.oicr.on.ca/pailab/private/xinghansun/UCSC_track_hub/humanFetalRhombicLip"


### Functions ###
## Utils ##
#' Prepares GRanges for export.bb() function. 
#' @details adds dummy score column and seqlengths() for provided genome
#' @param gr (GRanges) ranges to convert
#' @param genome (BSgenome) genome object. e.g., BSgenome.Hsapiens.UCSC.hg38.
prepareForBigBedExport <- function(gr,genome){
  gr$score <- 0
  gr <- sortSeqlevels(gr)
  tmp <- seqlengths(genome);
  tmp <- subset(tmp, names(tmp) %in% seqlevels(gr))
  tmp <- tmp[sortSeqlevels(names(tmp))]
  seqlengths(gr) <- tmp
  
  # theoretically, export.bb should take care of conversion btwn 1-based GRanges
  # and 0-based BigBed. However, it seems it doesn't. So have to make a pseudo
  # 0-based GRanges object for export.bb
  start(gr) <- start(gr) - 1
  
  return(gr)
}

#' Writes trackDb.txt entry for track.
#' @details Follows format from here: https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html#Setup
#' 
#' @param trackDbFile (char) relative path to trackDb.txt file to write to
#' @param hubUrl (char) URL to track hub directory
#' @param fileName (char) filename of track
#' @param trackName (char) track name. See UCSC doc.
#' @param shortLabel (char) shortLabel entry. See UCSC doc.
#' @param longLabel (char) longer descriptive name. See UCSC doc.
#' @param color (char) red,green,blue e.g., "255,0,0" for red. If NULL (default), is not included.
writeTrackDbEntry <- function(trackDbFile, hubUrl, trackName, fileName, shortLabel, longLabel,
                              color=NULL){
  message("writing track description")
  tFile <- trackDbFile
  cat(sprintf("track %s\n", trackName),
      file=tFile, append=TRUE)
  cat("type bigBed\n",file=tFile, append=TRUE)
  cat(sprintf("bigDataUrl %s/%s\n",hubUrl,fileName),
      file=tFile, append=TRUE)
  cat(sprintf("shortLabel %s\n",shortLabel),
      file=tFile, append=TRUE)
  cat(sprintf("longLabel %s\n", longLabel), 
      file=tFile, append=TRUE)
  if (!is.null(color)) {
    cat(sprintf("color %s\n", color),
        file=tFile,append=TRUE)
  }
  cat("\n",file=tFile,append=TRUE)
}


## tracks ##
# feCB H3K27ac peak summit #
addTrack_feCBsummit <- function(outDir, trackDbFile, srcUrl) {
  summit <- get_fetalCB_h3k27ac_summit(up = 1, down = 0)
  summit$Peak.Summit <- NULL
  
  summit <- prepareForBigBedExport(summit, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("feCB_H3K27ac_summit_hg38_%s.bb", dt)
  export.bb(summit, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="feCBk27acSummit_hg38",
    fileName=outFile,
    shortLabel="feCBk27acSummit_hg38",
    longLabel="Aldinger feCB H3K27ac peak summit",
    color="50,0,0"
  )
}


# dmr #
addTrack_dmr <- function(outDir, trackDbFile, srcUrl) {
  dmr <- get_cpg_dmrs()
  dmr_date <- get_configs("CPG_DMR_DATE")
  dmr$nCG <- NULL
  dmr$areaStat <- NULL
  dmr$name <- ifelse(dmr$diff.Methy>0, "RL-SVZ-hypo", "RL-SVZ-hyper")
  dmr$diff.Methy <- NULL
  
  dmr <- prepareForBigBedExport(dmr, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("rlDMR_hg38_%s.bb", dmr_date)
  export.bb(dmr, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="rlDMR_hg38",
    fileName=outFile,
    shortLabel="rlDMR_hg38",
    longLabel="DMRs between RL-SVZ and RL-VZ",
    color="50,0,0"
  )
}


# N2 #
addTrack_N2 <- function(outDir, trackDbFile, srcUrl) {
  N2_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                          "/ChIP-seq/Whalen_2023",
                          "/GSE110758_human-HS1-11-N2-pooled-K27ac.narrowPeak.gz")
  res <- rtracklayer::import(N2_h3k27ac_file)
  res <- liftOver_gr(res)

  res$pValue <- NULL
  res$qValue <- NULL
  res$peak <- NULL
  res$signalValue <- NULL
  res$score <- NULL
  res <- prepareForBigBedExport(res, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("N2_H3K27ac_hg38.bb")
  
  export.bb(res, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="N2_H3K27ac_hg38",
    fileName=outFile,
    shortLabel="N2_H3K27ac_hg38",
    longLabel="N2 H3K27ac pooled",
    color="50,0,0"
  )
}


# N3 #
addTrack_N3 <- function(outDir, trackDbFile, srcUrl) {
  N3_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                          "/ChIP-seq/Whalen_2023",
                          "/GSE110758_human-HS1-11-N3-pooled-K27ac.narrowPeak.gz")
  res <- rtracklayer::import(N3_h3k27ac_file)
  res <- liftOver_gr(res)
  
  res$pValue <- NULL
  res$qValue <- NULL
  res$peak <- NULL
  res$signalValue <- NULL
  res$score <- NULL
  res <- prepareForBigBedExport(res, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("N3_H3K27ac_hg38.bb")
  
  export.bb(res, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="N3_H3K27ac_hg38",
    fileName=outFile,
    shortLabel="N3_H3K27ac_hg38",
    longLabel="N3 H3K27ac pooled",
    color="50,0,0"
  )
}



# NPC K27ac #
addTrack_NPC27ac <- function(outDir, trackDbFile, srcUrl) {
  npc_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                           "/ChIP-seq/Choi_2021",
                           "/GSE158382_NPC_h3k27ac_merged_sorted_noDup_peaks.narrowPeak.gz")
  res <- rtracklayer::import(npc_h3k27ac_file)
  res <- liftOver_gr(res)
  
  res$pValue <- NULL
  res$qValue <- NULL
  res$peak <- NULL
  res$signalValue <- NULL
  res$score <- NULL
  res <- prepareForBigBedExport(res, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("NPC_H3K27ac_hg38.bb")
  
  export.bb(res, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="NPC_H3K27ac_hg38",
    fileName=outFile,
    shortLabel="NPC_H3K27ac_hg38",
    longLabel="NPC H3K27ac",
    color="50,0,0"
  )
}



# NPC ATAC #
addTrack_NPCatac <- function(outDir, trackDbFile, srcUrl) {
  npc_atac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                        "/ATAC-seq/Choi_2021",
                        "/GSE158382_NPC_ATAC_peaks.narrowPeak.gz")
  res <- rtracklayer::import(npc_atac_file)
  res <- liftOver_gr(res)

  res$pValue <- NULL
  res$qValue <- NULL
  res$peak <- NULL
  res$signalValue <- NULL
  res$score <- NULL
  
  res <- prepareForBigBedExport(res, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("NPC_ATAC_hg38.bb")
  
  export.bb(res, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="NPC_ATAC_hg38",
    fileName=outFile,
    shortLabel="NPC_ATAC_hg38",
    longLabel="NPC ATAC",
    color="50,0,0"
  )
}



# G34 SNV/INDEL #
addTrack_g34mut <- function(outDir, trackDbFile, srcUrl) {
  g34_mut <- import.bed(glue("/data/xsun/20240314/mutation/PEMECA-PCAWG/merged",
                             "/Group34_PEMECA-PCAWG_snv-indel_hg38.bed")
                        )
  
  g34_mut <- prepareForBigBedExport(g34_mut, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("g34_mut_hg38.bb")
  
  export.bb(g34_mut, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="g34_mut_hg38",
    fileName=outFile,
    shortLabel="g34_mut_hg38",
    longLabel="Group 3/4 MB SNV/INDEL",
    color="50,0,0"
  )
}


# G34 ATAC Smith et al #
addTrack_g34atac <- function(outDir, trackDbFile, srcUrl) {
  g3_atac <- import.bed("/data/xsun/20240314/cre/Smith2022-Group3/allEnhancerLikeElements.bed")
  g4_atac <- import.bed("/data/xsun/20240314/cre/Smith2022-Group4/allEnhancerLikeElements.bed")
  g3_atac$name <- "G3_ATAC"
  g4_atac$name <- "G4_ATAC"
  g34_atac <- c(g3_atac, g4_atac)
  # stick with common chrs
  g34_atac <- g34_atac[nchar(as.character(seqnames(g34_atac))) <= 5] 
  
  g34_atac <- prepareForBigBedExport(g34_atac, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("g34_atac_hg38.bb")
  
  export.bb(g34_atac, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="g34_atac_hg38",
    fileName=outFile,
    shortLabel="g34_atac_hg38",
    longLabel="Group 3/4 MB ATAC-seq combined peaks from Smith et al 2022",
    color="0,0,100"
  )
}



# ranked target #
addTrack_rankedTarget <- function(outDir, file, trackDbFile, srcUrl) {
  message(sprintf("getting ranked targets from: %s", file))
  
  targets <- read.table(file, stringsAsFactors = F, header = T, sep = "\t")
  targets_gr <- GenomicRanges::makeGRangesFromDataFrame(targets, ignore.strand = T)
  targets_gr$name <- paste0("rank-", targets$rank)

  targets_gr <- prepareForBigBedExport(targets_gr, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("rankedTargets_hg38_%s.bb", dt)
  export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="rankedTargets_hg38",
    fileName=outFile,
    shortLabel="rankedTargets_hg38",
    longLabel="RL DMRs MPRA ranked targets",
    color="50,0,0"
  )
}


addTrack_additionalTarget <- function(file, outDir, trackDbFile, srcUrl) {
  message(sprintf("getting additional targets from: %s", file))
  
  targets <- read.table(file, stringsAsFactors = F, header = T, sep = "\t")
  targets_gr <- GenomicRanges::makeGRangesFromDataFrame(targets, ignore.strand = T)

  targets_gr <- prepareForBigBedExport(targets_gr, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("additionalTargets_hg38_%s.bb", dt)
  export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="additionalTargets_hg38",
    fileName=outFile,
    shortLabel="additionalTargets_hg38",
    longLabel="MPRA additional targets",
    color="50,0,0"
  )
}


addTrack_posCtrl <- function(file, outDir, trackDbFile = tFile, srcUrl = srcUrl) {
  targets_gr <- import.bed(file)
  targets_gr$score <- NULL
  
  targets_gr <- prepareForBigBedExport(targets_gr, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("posCtrl_hg38_%s.bb", dt)
  export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="posCtrl_hg38",
    fileName=outFile,
    shortLabel="posCtrl_hg38",
    longLabel="MPRA ranked positive controls",
    color="200,0,0"
  )
}


addTrack_oligos <- function(file, outDir, trackDbFile = tFile, srcUrl = srcUrl) {
  targets_gr <- readRDS(file)
  
  targets_gr[targets_gr$class == "ranked"]$name <- 
    paste0("Rank-", targets_gr[targets_gr$class == "ranked"]$rank)
  
  targets_gr[targets_gr$class == "SNV"]$name <- 
    names(targets_gr[targets_gr$class == "SNV"])
  
  tmp <- targets_gr$name
  elementMetadata(targets_gr) <- NULL
  targets_gr$name <- tmp
  
  targets_gr <- prepareForBigBedExport(targets_gr, BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("designedOligos_hg38_%s.bb", dt)
  export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="designedOligos_hg38",
    fileName=outFile,
    shortLabel="designedOligos_hg38",
    longLabel="MPRA final oligos",
    color="200,0,0"
  )
}



addTrack_vistaHB <- function(outDir, trackDbFile = tFile, srcUrl = srcUrl) {
  ### VISTA hindbrain enhancer ###
  all_lines <- readLines(get_configs("VISTA_HB"))
  filtered_lines <- grep("^>", all_lines, value = TRUE)
  loc <- stringr::str_split(filtered_lines, "\\|", simplify = T)[,2]
  loc <- stringr::str_trim(loc)
  element <- stringr::str_split(filtered_lines, "\\|", simplify = T)[,3]
  element <- sub(" ", "_", stringr::str_trim(element))
  element <- glue("VISTA_{element}")
  
  vista_hb <- as.data.frame(stringr::str_split(loc, ":|-", simplify = T))
  colnames(vista_hb) <- c("seqnames", "start", "end")
  vista_hb_gr_hg19 <- GenomicRanges::makeGRangesFromDataFrame(vista_hb)
  vista_hb_gr_hg38 <- liftOver_gr(vista_hb_gr_hg19)
  vista_hb_gr_hg38$name <- element
  
  targets_gr <- prepareForBigBedExport(vista_hb_gr_hg38, 
                                       BSgenome.Hsapiens.UCSC.hg38)
  
  outFile <- sprintf("vistaHindbrain_hg38_%s.bb", dt)
  export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
  writeTrackDbEntry(
    trackDbFile=trackDbFile,
    hubUrl=srcUrl,
    trackName="vistaHindbrain_hg38",
    fileName=outFile,
    shortLabel="vistaHindbrain_hg38",
    longLabel="VISTA positive human hindbrain enhancers",
    color="0,50,0"
  )
}


#addTrack_g34sv_pcawg <- function(outDir, trackDbFile = tFile, srcUrl = srcUrl) {
#  # load hg38 sv
#  tmp <- read.table("/data/xsun/ol_test/wgs/pcawg_svs/hg38/Group34_sv_hg38.bed", 
#                    sep = "\t", header = F)
#  colnames(tmp) <- c("seqnames", "start", "end", "sample_id", "sv_class", "n_caller")
#  tmp$start <- tmp$start + 1 # bed to 1-based GRanges
#  pcawg_sv_hg38 <- GenomicRanges::makeGRangesFromDataFrame(tmp, keep.extra.columns = T)
#  
#  # add group info
#  id <- read.csv("/data/xsun/ol_test/wgs/pcawg_svs/id_subgroup.csv")
#  info <- id$Subgroup
#  names(info) <- id$id
#  pcawg_sv_hg38$subgroup <- info[pcawg_sv_hg38$sample_id]
#  
#  for (group in c("Group3", "Group4")) {
#    for (svclass in c("DEL", "DUP", "t2tINV", "h2hINV")) {
#      sub_gr <- pcawg_sv_hg38[pcawg_sv_hg38$subgroup == group & 
#                                pcawg_sv_hg38$sv_class == svclass
#                              ]
#      tmp <- split(sub_gr, sub_gr$sample_id)
#      sub_gr <- unlist(GRangesList(lapply(tmp, reduce)))
#      
#      sub_gr$name <- names(sub_gr)
#      names(sub_gr) <- NULL
#      
#      # add track
#      targets_gr <- prepareForBigBedExport(sub_gr, 
#                                           BSgenome.Hsapiens.UCSC.hg38)
#      
#      outFile <- sprintf("pcawg%s%s_hg38_%s.bb", group, svclass,dt)
#      export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
#      writeTrackDbEntry(
#        trackDbFile=trackDbFile,
#        hubUrl=srcUrl,
#        trackName= sprintf("pcawg%s%s_hg38", group, svclass),
#        fileName=outFile,
#        shortLabel= sprintf("pcawg%s%s_hg38", group, svclass),
#        longLabel=sprintf("PCAWG %s %s", group, svclass),
#        color="0,80,20"
#      )
#    }
#  }
#}

addTrack_g34sv_pcawg <- function(outDir, trackDbFile = tFile, srcUrl = srcUrl) {
  # load hg38 sv
  tmp <- read.table("/data/xsun/ol_test/wgs/pcawg_svs/hg38/Group34_sv_hg38.bed", 
                    sep = "\t", header = F)
  colnames(tmp) <- c("seqnames", "start", "end", "sample_id", "sv_class", "n_caller")
  tmp$start <- tmp$start + 1 # bed to 1-based GRanges
  pcawg_sv_hg38 <- GenomicRanges::makeGRangesFromDataFrame(tmp, keep.extra.columns = T)
  
  # add group info
  id <- read.csv("/data/xsun/ol_test/wgs/pcawg_svs/id_subgroup.csv")
  info <- id$Subgroup
  names(info) <- id$id
  pcawg_sv_hg38$subgroup <- info[pcawg_sv_hg38$sample_id]
  
  for (group in c("Group3", "Group4")) {
    for (svclass in c("DEL", "DUP", "t2tINV", "h2hINV")) {
      sub_gr <- pcawg_sv_hg38[pcawg_sv_hg38$subgroup == group & 
                                pcawg_sv_hg38$sv_class == svclass
      ]
      
      elementMetadata(sub_gr) <- elementMetadata(sub_gr)[,1]
      names(elementMetadata(sub_gr)) <- "name"
      
      # add track
      targets_gr <- prepareForBigBedExport(sub_gr, 
                                           BSgenome.Hsapiens.UCSC.hg38)
      
      outFile <- sprintf("pcawg%s%s_hg38_%s.bb", group, svclass,dt)
      export.bb(targets_gr, con = glue("{outDir}/{outFile}"))
      writeTrackDbEntry(
        trackDbFile=trackDbFile,
        hubUrl=srcUrl,
        trackName= sprintf("pcawg%s%s_hg38", group, svclass),
        fileName=outFile,
        shortLabel= sprintf("pcawg%s%s_hg38", group, svclass),
        longLabel=sprintf("PCAWG %s %s", group, svclass),
        color="0,80,20"
      )
    }
  }
}



### main ###
## peak summit ##
addTrack_feCBsummit(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)

## muts ##
addTrack_g34mut(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)

## dmr ##
addTrack_dmr(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)

## N2/N3 H3K27ac ##
addTrack_N2(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)
addTrack_N3(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)

## NPC K27ac/ATAC ##
addTrack_NPC27ac(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)
addTrack_NPCatac(outDir = outDir, trackDbFile = tFile, srcUrl = srcUrl)

## ranked targets ##
inFile <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/240711/oligoDesign",
               "/targetRanking_240731.tsv"
               )
addTrack_rankedTarget(
  file = inFile,
  outDir = outDir,
  trackDbFile = tFile,
  srcUrl = srcUrl
)

## additional targets ##
inFile <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/240711/oligoDesign",
               "/additionalTargets_240731.tsv"
)
addTrack_additionalTarget(
  file = inFile,
  outDir = outDir,
  trackDbFile = tFile,
  srcUrl = srcUrl
)


## pos control ##
inFile <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/240711/oligoDesign",
               "/posCtrlRanking_hg38_240802.bed")

addTrack_posCtrl(
  file = inFile,
  outDir = outDir,
  trackDbFile = tFile,
  srcUrl = srcUrl
)


## final oligos ##
inFile <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3",
               "/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711",
               "/oligoDesign/targetProbes_hg38_240731.rds")
addTrack_oligos(
  file = inFile,
  outDir = outDir,
  trackDbFile = tFile,
  srcUrl = srcUrl
)


## VISTA positive hindbrain enhancer ##
addTrack_vistaHB(
  outDir = outDir,
  trackDbFile = tFile,
  srcUrl = srcUrl
  )

## G34 ATAC from smith et al ##
addTrack_g34atac(
  outDir = outDir, 
  trackDbFile = tFile,
  srcUrl = srcUrl
)



## PCAWG G34 SV
addTrack_g34sv_pcawg(
  outDir = outDir,
  trackDbFile = tFile,
  srcUrl = srcUrl
)








