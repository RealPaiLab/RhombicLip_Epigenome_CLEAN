rm(list=ls())

library(Biostrings)
library(universalmotif)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggpubr)

setwd(this.path::this.dir())
source("./utils.R")
source("../../EM-seq/Analysis/FETHB2-FETHB3/utils.R")

options(
  meme_bin = "/opt/bin"
)

dt <- format(Sys.Date(),"%y%m%d")
dmr_date <- get_configs("CPG_DMR_DATE")
outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/{dmr_date}/oligoDesign")

configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN

whalenNegCtrl_file <- sprintf("%s/whalen_MPRA_selected_negCtrl_hg19_%s.bed", 
                          outDir, "240802")

logFile <- sprintf("%s/negCtrl_ranking_%s.log",
                   outDir, 
                   dt
)

# use targets to generate shuffled neg ctrls
bedFile <- sprintf("%s/targetProbes_hg38_%s.bed", outDir, "240731") 

configs <- yaml::read_yaml("./config.yaml")
budget <- configs$N_NEG


main <- function() {
  message(sprintf("--- Loading %s as source of shuffling.", bedFile))
  bed <- import.bed(bedFile)
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  seq <- getSeq(genome, bed)
  names(seq) <- paste0("seq_", 1:length(seq))
  name_dict <- bed$name
  names(name_dict) <- paste0("seq_", 1:length(seq))
  
  message("--- Shuflling sequences")
  shuffled_seq <- shuffle_sequences(seq, 
                                    k = 1, 
                                    #method = "euler", 
                                    nthreads = 4,
                                    rng.seed = 1234
  )
  
  message("--- Running MEME FIMO on the unshuffled sequences")
  original_fimo <- memes::runFimo(seq, get_configs("HOCOMOCO_V12"))
  
  message("--- Running MEME FIMO on the shuffled sequences")
  shuffled_fimo <- memes::runFimo(shuffled_seq, get_configs("HOCOMOCO_V12"))
  
  o_df <- as.data.frame(table(original_fimo@seqnames)) %>% 
    mutate(sequence = "Raw")
  s_df <- as.data.frame(table(shuffled_fimo@seqnames)) %>%
    mutate(sequence = "Shuffled")
  os_df <- rbind(o_df, s_df)
  
  message("--- Saving FIMO motif count boxplot comparing shuffled and unshuffled.")
  .plot <- ggplot(os_df, aes(x = sequence, y = Freq, fill = sequence)) + 
    geom_violin() + 
    scale_fill_manual(values = c(Raw = "royalblue", Shuffled = "#FF7F50")) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    stat_compare_means(method = "wilcox.test")
  outFile <- sprintf("%s/shuffledSeq_FIMO_violinPlot_%s.png", outDir, dt)
  ggsave(outFile, .plot, dpi = 600, height = 8, width = 10, units = "in")
  
  
  # select shuffled sequences with lowest number of FIMO motifs
  message("Ranking shuffled sequences based on number of FIMO motifs")
  ranked_intervals <- s_df %>% 
    arrange(Freq) %>% 
    pull(Var1)
  shuffled_seq_ranked <- shuffled_seq[ranked_intervals]
  names(shuffled_seq_ranked) <- name_dict[names(shuffled_seq_ranked)]
  names(shuffled_seq_ranked) <- paste0("shuffled_", names(shuffled_seq_ranked))
  
  ### add Whalen 2023 selected neg controls
  negCtrl_gr <- import.bed(whalenNegCtrl_file)
  negCtrl_gr <- liftOver_gr(negCtrl_gr)
  resized_negCtrl_gr <- resize_interval(negCtrl_gr, mpra_target_len)
  
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  resized_negCtrl_seqs <- getSeq(genome, resized_negCtrl_gr)
  names(resized_negCtrl_seqs) <- negCtrl_gr$name

  # bind as top ranked
  final <- c(resized_negCtrl_seqs, shuffled_seq_ranked)

  ### output top seqs ###
  shuffled_seq_ranked <- final[1:budget]
  
  check_SceI(shuffled_seq_ranked)
  
  outFile <- sprintf("%s/negCtrlRanking_hg38_%s.fasta", outDir, dt)
  message(sprintf("Writing sequences to %s", outFile))
  memes::write_fasta(shuffled_seq_ranked, outFile)
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


