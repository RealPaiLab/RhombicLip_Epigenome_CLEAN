# run TFBS motif enrichment test using AME on DMRs identified between RL-SVZ and RL-VZ

rm(list=ls())

library(memes) # MEME v5.5.5
library(BSgenome.Hsapiens.UCSC.hg38)
library(glue)
source("../../utils.R")

dt <- format(Sys.Date(),"%y%m%d")
tf_db <- get_configs("HOCOMOCO_V12")

#"/home/rstudio/software/MEME/motif_databases/HUMAN/HOCOMOCOv11_CORE_HUMAN_mono_meme_format.meme"
meme_path <- "/home/rstudio/software/meme/bin"

rootDir <- glue("/home/rstudio/isilon/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")
outDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection/%s",
                  rootDir, 
                  get_configs("CPG_DMR_DATE")
                  )
dmrFile <- sprintf("%s/moreStringent_251210/DMRs_moreStringent.csv", outDir)                

outDir <- dirname(dmrFile)

logFile <- sprintf("%s/AME_fullMEME_%s.log",outDir,dt)
sink(logFile, split=TRUE)

tryCatch({
  message("Motif enrichment test on DMRs using AME")
  message(sprintf("Output will be written under: %s\n",outDir))
  
  ### Load DMRs ###
  ranges_dmrs <- read.delim(dmrFile,header=TRUE,sep=",",
                            stringsAsFactors = FALSE
                            )
  
  ### All dmr ###
  all_dmr <- ranges_dmrs
    all_dmr <- GRanges(
        seqnames = all_dmr$chr,
        ranges = IRanges(start=all_dmr$start,
                         end=all_dmr$end
                         )
        )
  message(sprintf("%s DMRs to be analyzed",length(all_dmr)))
  all_dmr <- GenomicRanges::reduce(all_dmr)
  all_dmr_seq <- get_sequence(all_dmr, BSgenome.Hsapiens.UCSC.hg38)
  
  all_dmr_outDir <- sprintf("%s/all_dmr_AME_fullMEME_%s", outDir, dt)
  
 options(meme_bin = meme_path)
 browser()
  all_dmr_ame <- runAme(
    input = all_dmr_seq, 
    outdir = all_dmr_outDir, 
    database = tf_db
    )
  message(sprintf("%i TFs identified", nrow(all_dmr_ame)))
  message(sprintf("Top TFs enriched in the regions are: \n- %s", 
                  paste(head(all_dmr_ame$motif_id, 10), 
                        collapse = "\n- "
                        )
                  )
          )

    x <- as.data.frame(all_dmr_ame)
}, error = function(e) {
    message("Error encountered:")
    print(e)
}, finally = {
    sink()
})


