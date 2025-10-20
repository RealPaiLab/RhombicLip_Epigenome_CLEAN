rm(list=ls())

library(GenomicRanges)
library(rtracklayer)
library(glue)


default_wd <- this.path::this.dir()
setwd(default_wd) # set current scripts' dir as working dir
source("../../EM-seq/Analysis/FETHB2-FETHB3/utils.R")

dt <- format(Sys.Date(),"%y%m%d")
dmr_date <- get_configs("CPG_DMR_DATE")
outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/{dmr_date}/oligoDesign")

logFile <- sprintf("%s/generate_finalOligoPool_%s.log",
                   outDir, 
                   dt
)

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN

# fasta file paths
crs_fasta <- sprintf("%s/targetProbes_hg38_%s.fasta", outDir, 240731)
pos_fasta <- sprintf("%s/posCtrlRanking_hg38_%s.fasta", outDir, 240802)
neg_fasta <- sprintf("%s/negCtrlRanking_hg38_%s.fasta", outDir, 240802)

#' add adapters to the provided sequences
#' @param seq (DNAStringSet) target sequences
#' @return (DNAStringSet)
add_adapter <- function(seq) {
  # availabel adpaters from Gordon 2020
  # 5′-AGGACCGGATCAACT**200base_CRS**CATTGCGTGAACCGA-3′
  # 5′-AATGCTAGCGCATGG**200base_CRS**CTGCAACCTACGGAA-3′
  # 5′-TTACGAGCCGTAGTC**200base_CRS**GCATCTCAACGTGGT-3′
  
  adapterPair <- c("AGGACCGGATCAACT", "CATTGCGTGAACCGA")
  
  
  res <- DNAStringSet(paste0(adapterPair[1], as.character(seq), adapterPair[2]))
  names(res) <- names(seq)
  
  return(res)
}


#' generate 3-column format for Agilent SureDesign
#' @param seqs (DNAStringSet)
#' @return (data.frame)
generate_3col_agilent <- function(seqs) {
  # three columns SequenceName	Sequence	Replication	
  # SequenceName must have a len <= 100
  res <- data.frame(SequenceName = substr(names(seqs), 1, 100),
                    Sequence = as.character(seqs),
                    Replication = 1
                    )
  
  return(res)
}


main <- function() {
  crs <- readDNAStringSet(crs_fasta)
  pos <- readDNAStringSet(pos_fasta)
  neg <- readDNAStringSet(neg_fasta)
  
  names(pos) <- paste0("POSITIVE_", names(pos))
  names(neg) <- paste0("NEGATIVE_", names(neg))
  
  pool <- c(crs, pos, neg)
  pool_with_adapter <- add_adapter(pool)
  
  message(sprintf("Total %d sequences (CRS: %d; POS: %d; NEG: %d)", 
                  length(pool_with_adapter), 
                  length(crs), length(pos), length(neg)
                  )
          )
  
  message(sprintf("Median test region size is %d (min: %d; max: %d)",
                  median(width(pool)),
                  min(width(pool)), 
                  max(width(pool))
  ))
  
  message(sprintf("Median oligo size is %d (min: %d; max: %d)",
                  median(width(pool_with_adapter)),
                  min(width(pool_with_adapter)), 
                  max(width(pool_with_adapter))
                  ))
  
  # fasta format
  outFile <- sprintf("%s/finalOligos_hg38_%s.fasta", outDir, dt)
  message(sprintf("Writing final oligos in fasta format to %s", outFile))
  memes::write_fasta(pool_with_adapter, outFile)
  
  # agilent 3-col format
  # The file extension must be .list or .csv or .txt or .tdt.
  outFile <- sprintf("%s/finalOligos_hg38_%s.list", outDir, dt)
  message(sprintf("Writing final oligos in 3-col Agilent format to %s", outFile))
  ag <- generate_3col_agilent(pool_with_adapter)
  write.table(ag, outFile, col.names = T, row.names = F, quote = F, sep = "\t")
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
