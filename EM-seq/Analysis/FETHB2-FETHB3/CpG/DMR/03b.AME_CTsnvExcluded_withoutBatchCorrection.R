rm(list=ls())

library(memes) # MEME v5.5.5
library(BSgenome.Hsapiens.UCSC.hg38)
library(glue)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")

### config ###

options(
  meme_bin = "/opt/bin"
)

dt <- format(Sys.Date(),"%y%m%d")
tf_db <- get_configs("TF_DB")

rootDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")
outDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection/%s",
                  rootDir, 
                  get_configs("CPG_DMR_DATE")
                  )
logFile <- sprintf("%s/AME_CTsnvExcluded_withoutBatchCorrection_%s_%s.log",
                   outDir, 
                   tf_db, 
                   dt)


main <- function() {
  message("Motif enrichment test on DMRs using AME")
  message(sprintf("Output will be written under: %s\n",outDir))
  
  ### Get RL HOCOMOCOv12 db for meme ###
  rl_meme_db <- get_configs("RL_MEME_DB")
  
  ### Load DMRs ###
  ranges_dmrs <- get_cpg_dmrs()
  
  ### RL-SVZ diff hypo ###
  if (sum(ranges_dmrs$diff.Methy > 0) > 0) {
    svz_diff_hypo <- subset(ranges_dmrs, diff.Methy > 0)
    message(sprintf("%s DMRs are hypomethylated in RL-SVZ comparing to RL-VZ", 
                    length(svz_diff_hypo)
                    )
            )
    svz_diff_hypo <- GenomicRanges::reduce(svz_diff_hypo)
    svz_diff_hypo_seq <- get_sequence(svz_diff_hypo, 
                                      BSgenome.Hsapiens.UCSC.hg38
                                      )
    
    ## only normal RL TFs ##
    svz_diff_hypo_outDir <- glue("{outDir}/SVZ_diff_hypo_AME_activeInRL_{tf_db}_{dt}")

    message(glue("Running AME, results will be saved under: ",
                 "{svz_diff_hypo_outDir}"
                 )
            )
    message(sprintf("Using meme db: %s", rl_meme_db))
    options(meme_db = rl_meme_db)
    
    svz_diff_hypo_ame <- runAme(input = svz_diff_hypo_seq, 
                                outdir = svz_diff_hypo_outDir
                                )
    message(sprintf("%i TFs identified", nrow(svz_diff_hypo_ame)))
    message(
      sprintf(
        "Top TFs enriched in the regions are: \n- %s",
        paste(head(svz_diff_hypo_ame$motif_id, 10), collapse = "\n- ")
        )
      )
  } else {
    message("No DMRs are hypomethylated in RL-SVZ comparing to RL-VZ.")
  }
  
  ### RL-SVZ diff hyper ###
  if (sum(ranges_dmrs$diff.Methy < 0) > 0) {
    svz_diff_hyper <- subset(ranges_dmrs, diff.Methy < 0)
    message(sprintf("%s DMRs are hypermethylated in RL-SVZ comparing to RL-VZ",
                    length(svz_diff_hyper)
                    )
            )
    svz_diff_hyper <- GenomicRanges::reduce(svz_diff_hyper)
    svz_diff_hyper_seq <- get_sequence(svz_diff_hyper, 
                                       BSgenome.Hsapiens.UCSC.hg38
                                       )
    
    ## only normal RL TFs ##
    svz_diff_hyper_outDir <- sprintf("%s/SVZ_diff_hyper_AME_activeInRL_%s_%s", 
                                     outDir,
                                     tf_db,
                                     dt)
    
    message(glue("Running AME, results will be saved under: ",
                 "{svz_diff_hyper_outDir}"
                 )
            )
    message(sprintf("Using meme db: %s", rl_meme_db))
    options(meme_db = rl_meme_db)
    
    svz_diff_hyper_ame <- runAme(input = svz_diff_hyper_seq, 
                                 outdir = svz_diff_hyper_outDir
                                 )
    message(sprintf("%i TFs identified", nrow(svz_diff_hyper_ame)))
    message(sprintf("Top TFs enriched in the regions are: \n- %s",
                    paste(
                      head(svz_diff_hyper_ame$motif_id, 10), 
                      collapse = "\n- "
                      )
                    )
            )
  } else {
    message("No DMRs are hypermethylated in RL-SVZ comparing to RL-VZ.")
  }
  
  ### All dmr ###
  all_dmr <- ranges_dmrs
  message(sprintf("%s DMRs to be analyzed",length(all_dmr)))
  all_dmr <- GenomicRanges::reduce(all_dmr)
  all_dmr_seq <- get_sequence(all_dmr, BSgenome.Hsapiens.UCSC.hg38)
  
  ## only normal RL TFs ##
  all_dmr_outDir <- glue("{outDir}/all_dmr_AME_activeInRL_{tf_db}_{dt}")
  
  message(glue("Running AME, results will be saved under: ",
               "{all_dmr_outDir}"
               )
          )
  message(sprintf("Using meme db: %s", rl_meme_db))
  options(meme_db = rl_meme_db)
  
  all_dmr_ame <- runAme(input = all_dmr_seq, outdir = all_dmr_outDir)
  message(sprintf("%i TFs identified", nrow(all_dmr_ame)))
  message(sprintf("Top TFs enriched in the regions are: \n- %s", 
                  paste(head(all_dmr_ame$motif_id, 10), 
                        collapse = "\n- "
                        )
                  )
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
         warning = function(w) {message(w)},
         finally = {
           message("\n\n--------- R sessionInfo ---------\n\n")
           print(sessionInfo())
           sink(type = "output")
           sink(type = "message")
           close(logFileCon)
         }
) 

























