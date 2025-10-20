rm(list=ls())

library(memes) # MEME v5.5.5
library(BSgenome.Hsapiens.UCSC.hg38)

### config ###

options(meme_db = "/data/xsun/db/meme/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme",
        meme_bin = "/opt/bin"
        )

dt <- format(Sys.Date(),"%y%m%d")

rootDir <- "/data/xsun/output/EMseq_FETHB3"
outDir <- sprintf("%s/DMRs/CTsnv_included/%s",rootDir, "240330")

dmrFile <- "/data/xsun/output/EMseq_FETHB3/DMRs/CTsnv_included/240330/DMRs.csv"

logFile <- sprintf("%s/AME_%s.log",outDir, dt)


main <- function() {
  df_dmrs <- read.table(dmrFile, stringsAsFactors = F, header = T)
  ranges_dmrs <- GRanges(df_dmrs$chr, IRanges(df_dmrs$start, df_dmrs$end))
  ranges_dmrs$diff.Methy <- df_dmrs$diff.Methy

  ### RL-SVZ diff hypo
  svz_diff_hypo <- subset(ranges_dmrs, diff.Methy > 0)
  message(sprintf("%s DMRs are hypomethylated in RL-SVZ comparing to RL-VZ", length(svz_diff_hypo)))
  svz_diff_hypo <- GenomicRanges::reduce(svz_diff_hypo)
  svz_diff_hypo_seq <- get_sequence(svz_diff_hypo, BSgenome.Hsapiens.UCSC.hg38)
  svz_diff_hypo_outDir <- sprintf("%s/SVZ_diff_hypo_AME_%s", outDir, dt)
  print(svz_diff_hypo_outDir)
  message(sprintf("Running AME, results will be saved under: %s", svz_diff_hypo_outDir))
  
  svz_diff_hypo_ame <- runAme(input = svz_diff_hypo_seq, outdir = svz_diff_hypo_outDir)
  message(sprintf("%i TFs identified", nrow(svz_diff_hypo_ame)))
  message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(svz_diff_hypo_ame$motif_id, 10), collapse = "\n- ")))
  
  ### RL-SVZ diff hyper
  svz_diff_hyper <- subset(ranges_dmrs, diff.Methy < 0)
  message(sprintf("%s DMRs are hypermethylated in RL-SVZ comparing to RL-VZ",length(svz_diff_hyper)))
  svz_diff_hyper <- GenomicRanges::reduce(svz_diff_hyper)
  svz_diff_hyper_seq <- get_sequence(svz_diff_hyper, BSgenome.Hsapiens.UCSC.hg38)
  svz_diff_hyper_outDir <- sprintf("%s/SVZ_diff_hyper_AME_%s", outDir, dt)
  message(sprintf("Running AME, results will be saved under: %s", svz_diff_hyper_outDir))
  
  svz_diff_hyper_ame <- runAme(input = svz_diff_hyper_seq, outdir = svz_diff_hyper_outDir)
  message(sprintf("%i TFs identified", nrow(svz_diff_hyper_ame)))
  message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(svz_diff_hyper_ame$motif_id, 10), collapse = "\n- ")))
  
  ### All dmr
  all_dmr <- ranges_dmrs
  message(sprintf("%s DMRs to be analyzed",length(all_dmr)))
  all_dmr <- GenomicRanges::reduce(all_dmr)
  all_dmr_seq <- get_sequence(all_dmr, BSgenome.Hsapiens.UCSC.hg38)
  all_dmr_outDir <- sprintf("%s/all_dmr_AME_%s", outDir, dt)
  message(sprintf("Running AME, results will be saved under: %s", all_dmr_outDir))
  
  all_dmr_ame <- runAme(input = all_dmr_seq, outdir = all_dmr_outDir)
  message(sprintf("%i TFs identified", nrow(all_dmr_ame)))
  message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(all_dmr_ame$motif_id, 10), collapse = "\n- ")))
}



### main ###
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

























