rm(list=ls())

library(memes) # MEME v5.5.5
library(BSgenome.Hsapiens.UCSC.hg38)

### config ###

options(
  meme_bin = "/opt/bin"
)

dt <- format(Sys.Date(),"%y%m%d")

# vz_meme_db <- "/.mounts/labs/pailab/private/xsun/output/RL_genes/20240423/filteredTFdb/Hendrikse2022_RLVZ_activeGenes_TF.meme"
# svz_meme_db <- "/.mounts/labs/pailab/private/xsun/output/RL_genes/20240423/filteredTFdb/Hendrikse2022_RLSVZ_activeGenes_TF.meme"
rl_meme_db <- "/.mounts/labs/pailab/private/xsun/output/RL_genes/20240423/filteredTFdb/Hendrikse2022_RL_activeGenes_TF.meme"

# vz_g34_meme_db <- "/.mounts/labs/pailab/private/xsun/output/RL_genes/20240423/filteredTFdb/RLVZ_G34_activeGenes_TF.meme"
# svz_g34_meme_db <- "/.mounts/labs/pailab/private/xsun/output/RL_genes/20240423/filteredTFdb/RLSVZ_G34_activeGenes_TF.meme"
rl_g34_meme_db <- "/.mounts/labs/pailab/private/xsun/output/RL_genes/20240423/filteredTFdb/RL_G34_activeGenes_TF.meme"


rootDir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG"
outDir <- sprintf("%s/DMRs/CTsnv_included/withoutBatchCorrection/%s",rootDir, "240419")

dmrFile <- sprintf("%s/DMRs.csv", outDir )

logFile <- sprintf("%s/AME_dssBatchCorrection_%s.log",outDir, dt)


main <- function() {
  df_dmrs <- read.table(dmrFile, stringsAsFactors = F, header = T)
  ranges_dmrs <- GRanges(df_dmrs$chr, IRanges(df_dmrs$start, df_dmrs$end))
  ranges_dmrs$areaStat <- df_dmrs$areaStat
  
  ### RL-SVZ diff hypo
  if (sum(ranges_dmrs$areaStat > 0) > 0) {
    svz_diff_hypo <- subset(ranges_dmrs, areaStat > 0)
    message(sprintf("%s DMRs are hypomethylated in RL-SVZ comparing to RL-VZ", length(svz_diff_hypo)))
    svz_diff_hypo <- GenomicRanges::reduce(svz_diff_hypo)
    svz_diff_hypo_seq <- get_sequence(svz_diff_hypo, BSgenome.Hsapiens.UCSC.hg38)
    
    ## only normal RL
    svz_diff_hypo_outDir <- sprintf("%s/SVZ_diff_hypo_AME_activeInRL_%s", outDir, dt)

    message(sprintf("Running AME, results will be saved under: %s", svz_diff_hypo_outDir))
    message(sprintf("Using meme db: %s", rl_meme_db))
    options(meme_db = rl_meme_db)
    
    svz_diff_hypo_ame <- runAme(input = svz_diff_hypo_seq, outdir = svz_diff_hypo_outDir)
    message(sprintf("%i TFs identified", nrow(svz_diff_hypo_ame)))
    message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(svz_diff_hypo_ame$motif_id, 10), collapse = "\n- ")))
    
    ## RL+tumour
    svz_diff_hypo_outDir <- sprintf("%s/SVZ_diff_hypo_AME_activeInRLg34_%s", outDir, dt)
    
    message(sprintf("Running AME, results will be saved under: %s", svz_diff_hypo_outDir))
    message(sprintf("Using meme db: %s", rl_g34_meme_db))
    options(meme_db = rl_g34_meme_db)
    
    svz_diff_hypo_ame <- runAme(input = svz_diff_hypo_seq, outdir = svz_diff_hypo_outDir)
    message(sprintf("%i TFs identified", nrow(svz_diff_hypo_ame)))
    message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(svz_diff_hypo_ame$motif_id, 10), collapse = "\n- ")))
  } else {
    message("No DMRs are hypomethylated in RL-SVZ comparing to RL-VZ.")
  }
  
  ### RL-SVZ diff hyper
  if (sum(ranges_dmrs$areaStat < 0) > 0) {
    svz_diff_hyper <- subset(ranges_dmrs, areaStat < 0)
    message(sprintf("%s DMRs are hypermethylated in RL-SVZ comparing to RL-VZ",length(svz_diff_hyper)))
    svz_diff_hyper <- GenomicRanges::reduce(svz_diff_hyper)
    svz_diff_hyper_seq <- get_sequence(svz_diff_hyper, BSgenome.Hsapiens.UCSC.hg38)
    
    ## only normal RL
    svz_diff_hyper_outDir <- sprintf("%s/SVZ_diff_hyper_AME_activeInRL_%s", outDir, dt)
    
    message(sprintf("Running AME, results will be saved under: %s", svz_diff_hyper_outDir))
    message(sprintf("Using meme db: %s", rl_meme_db))
    options(meme_db = rl_meme_db)
    
    svz_diff_hyper_ame <- runAme(input = svz_diff_hyper_seq, outdir = svz_diff_hyper_outDir)
    message(sprintf("%i TFs identified", nrow(svz_diff_hyper_ame)))
    message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(svz_diff_hyper_ame$motif_id, 10), collapse = "\n- ")))
    
    ## RL+tumour
    svz_diff_hyper_outDir <- sprintf("%s/SVZ_diff_hyper_AME_activeInRLg34_%s", outDir, dt)
    
    message(sprintf("Running AME, results will be saved under: %s", svz_diff_hyper_outDir))
    message(sprintf("Using meme db: %s", rl_g34_meme_db))
    options(meme_db = rl_g34_meme_db)
    
    svz_diff_hyper_ame <- runAme(input = svz_diff_hyper_seq, outdir = svz_diff_hyper_outDir)
    message(sprintf("%i TFs identified", nrow(svz_diff_hyper_ame)))
    message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(svz_diff_hyper_ame$motif_id, 10), collapse = "\n- ")))
  } else {
    message("No DMRs are hypermethylated in RL-SVZ comparing to RL-VZ.")
  }
  
  ### All dmr
  all_dmr <- ranges_dmrs
  message(sprintf("%s DMRs to be analyzed",length(all_dmr)))
  all_dmr <- GenomicRanges::reduce(all_dmr)
  all_dmr_seq <- get_sequence(all_dmr, BSgenome.Hsapiens.UCSC.hg38)
  
  ## RL only
  all_dmr_outDir <- sprintf("%s/all_dmr_AME_activeInRL_%s", outDir, dt)
  
  message(sprintf("Running AME, results will be saved under: %s", all_dmr_outDir))
  message(sprintf("Using meme db: %s", rl_meme_db))
  options(meme_db = rl_meme_db)
  
  all_dmr_ame <- runAme(input = all_dmr_seq, outdir = all_dmr_outDir)
  message(sprintf("%i TFs identified", nrow(all_dmr_ame)))
  message(sprintf("Top TFs enriched in the regions are: \n- %s", paste(head(all_dmr_ame$motif_id, 10), collapse = "\n- ")))
  
  ## RL+tumour
  all_dmr_outDir <- sprintf("%s/all_dmr_AME_activeInRLg34_%s", outDir, dt)
  
  message(sprintf("Running AME, results will be saved under: %s", all_dmr_outDir))
  message(sprintf("Using meme db: %s", rl_g34_meme_db))
  options(meme_db = rl_g34_meme_db)
  
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

























