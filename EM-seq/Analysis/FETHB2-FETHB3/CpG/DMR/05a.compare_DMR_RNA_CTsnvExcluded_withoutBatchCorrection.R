### NOTE: if you encounter problem saving the pdf figures, it's most likely 
### because of the warnings from ggrepel, which is normal, but will cause 
### tryCatch to stop plotting. Please run scripts in main function manually.
rm(list=ls())

library(glue)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)


setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")
source("../../CpG/overlapEnrichment/utils.R")
source("../../../../../MPRA/oligo_design/utils.R")

dt <- format(Sys.Date(),"%y%m%d")
cpg_dmr_date <- get_configs("CPG_DMR_DATE")

rootDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")

outDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection/%s/DMR_vs_RNA",
                  rootDir, 
                  cpg_dmr_date
                  )

logFile <- sprintf("%s/compare_DMR_RNA_CTsnvExcluded_withoutBatchCorrection_%s.log",
                   outDir, 
                   dt
)



plot_meth_rna_boxplot <- function(df, title = "") {
  .plot <- ggplot(df, aes(x = upRegIn, y = diff.Methy*100, fill = upRegIn)) + 
    geom_boxplot() +
    scale_fill_manual(name = "", 
                      values = c("RL-VZ" = "#fec140", "RL-SVZ" = "#88b728")
                      ) +
    geom_jitter(position = position_jitter(seed = 1234), alpha = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", alpha = 0.5) +
    stat_compare_means(method = "wilcox.test", label.x = 1, size = 8) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 22), 
      axis.text = element_text(size = 25),  
      legend.text = element_text(size = 22),
      legend.key.size = unit(1.5, "cm")
    ) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, 
             alpha = 0.05, fill = "lightcoral"
             ) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, 
             alpha = 0.05, fill = "lightblue"
             ) +
    ggtitle(title) +
    ylab("DMR methylation VZ - SVZ (%)") +
    xlab("") +
    scale_x_discrete(labels = c(
      "RL-VZ" = sprintf("Overexpressed in\nRL-VZ\n(N = %d)", 
                        sum(df$upRegIn == "RL-VZ")
                        ), 
      "RL-SVZ" = sprintf("Overexpressed in\nRL-SVZ\n(N = %d)", 
                         sum(df$upRegIn == "RL-SVZ")
                         )
      ))
    
  return(.plot)
}


plot_meth_rna_scatter <- function(df, 
                                  pearson_cor = TRUE,
                                  label_col = "nearestGene",
                                  logFC_thresh = 1.0,
                                  diff.Methy_thresh = 0.1,
                                  title = ""
                                  ) {
  scale <- ifelse(pearson_cor, 1.4, 1.15)
  .plot <- ggplot(df, aes(x = logFC, y = diff.Methy, label = nearestGene)) + 
    # highlight points
    geom_point(
      color = case_when(
        df$logFC <= -1*logFC_thresh & df$diff.Methy >= diff.Methy_thresh ~ "#88b728",
        df$logFC >= logFC_thresh & df$diff.Methy <= -1*diff.Methy_thresh ~ "#fec140",
        df$logFC <= -1*logFC_thresh & df$diff.Methy <= -1*diff.Methy_thresh ~ "royalblue",
        df$logFC >= logFC_thresh & df$diff.Methy >= diff.Methy_thresh ~ "royalblue",
        .default = "grey50"
        ),
      size = ifelse(
        abs(df$logFC) >= logFC_thresh & abs(df$diff.Methy) >= diff.Methy_thresh,
        1,
        0.2
        )
      ) +
    # label points
    geom_text_repel(
      max.overlaps  = 15,
      force = 2,
      data = subset(
        df, 
        abs(logFC) >= logFC_thresh & abs(diff.Methy) >= diff.Methy_thresh
        ),
      min.segment.length = 0, 
      segment.size = 0.2,
      segment.color = "grey50", 
      size = 4,
      seed = 1234
    ) +
    theme_classic() +
    # mark thresholds
    geom_hline(yintercept = diff.Methy_thresh, 
               linetype = "dashed", color = "grey", alpha = 0.5) +
    geom_hline(yintercept = -1*diff.Methy_thresh, 
               linetype = "dashed", color = "grey", alpha = 0.5) +
    geom_vline(xintercept = logFC_thresh, 
               linetype = "dashed", color = "grey", alpha = 0.5) +
    geom_vline(xintercept = -1*logFC_thresh, 
               linetype = "dashed", color = "grey", alpha = 0.5) +
    # annotate RNA tissue types
    annotate("segment", 
             x = min(df$logFC)+1, xend = min(df$logFC), 
             y = max(df$diff.Methy)*scale, yend = max(df$diff.Methy)*scale, 
             arrow = arrow(type = "closed"), color = "#88b728", linewidth = 1,) +
    annotate("segment", 
             x = max(df$logFC)-1, xend = max(df$logFC), 
             y = max(df$diff.Methy)*scale, yend = max(df$diff.Methy)*scale, 
             arrow = arrow(type = "closed"), color = "#fec140", linewidth = 1) +
    # annotate DNAm tissue types
    annotate("segment", 
             x = min(df$logFC)*1.15, xend = min(df$logFC)*1.15, 
             y = max(df$diff.Methy)*scale-0.1, yend = max(df$diff.Methy)*scale, 
             arrow = arrow(type = "closed"), color = "#88b728", linewidth = 1,) +
    annotate("segment", 
             x = min(df$logFC)*1.15, xend = min(df$logFC)*1.15, 
             y = min(df$diff.Methy)+0.1, yend = min(df$diff.Methy), 
             arrow = arrow(type = "closed"), color = "#fec140", linewidth = 1) +
    ggtitle(title) +
    ylab("DMR methylation VZ - SVZ (%)") + 
    xlab("RNA log2 fold-change VZ/SVZ")
  
  if (pearson_cor) {
    .plot <- .plot + 
      geom_smooth(method = "lm", col = alpha("red", 0.2), alpha = 0.1) +
      stat_cor(method = "pearson", 
               label.x = min(df$logFC), 
               label.y = max(df$diff.Methy)*1.2
               )
  }
  
  return(.plot)
}


main <- function() {
  ### Load data ###
  ## RL-VZ vs SVZ LCM DMR delta ##
  # Methyl1 = VZ; Methyl2 = SVZ; diff = VZ - SVZ
  dmr <- read.table(get_configs("CPG_DMR_FILE"), 
                    stringsAsFactors = F, 
                    header = T
                    )
  dmr_gr <- GenomicRanges::makeGRangesFromDataFrame(dmr, keep.extra.columns = T)
  
  ## RL-VZ vs SVZ LCM RNA logFC ##
  # logFC = VZ/SVZ
  deg <- read.table(get_configs("RL_DEG"), stringsAsFactors = F, header = T)
  deg$geneName <- rownames(deg)
  deg$upRegIn <- ifelse(deg$logFC > 0, "RL-VZ", "RL-SVZ")

  ### Link dmr and deg ###
  message("Annotating DMRs with nearest gene")
  #dmr_gr_annotated <- getNearestGene(dmr_gr)
  
  message("Separating DMRs by overlapping with promoters.")
  promoters_gr <- get_promoters() # default TSS +/- 1kb 
  ol <- findOverlaps(dmr_gr, promoters_gr)
  promoterDMR_gr <- dmr_gr[queryHits(ol)]
  promoterDMR_gr$nearestGene <- promoters_gr[subjectHits(ol)]$name
  promoterDMR_gr <- unique(promoterDMR_gr)
  
  nonpromoterDMR_gr <- getNearestGene(dmr_gr[-queryHits(ol)])
  message(sprintf("- Promoter DMRs: %i; non-promoter DMRs: %i",
                  length(promoterDMR_gr),
                  length(nonpromoterDMR_gr)
                  )
          )
  
  promoterDMR_df <- as.data.frame(promoterDMR_gr)
  nonpromoterDMR_df <- as.data.frame(nonpromoterDMR_gr)
  
  promoterDMR_df <- merge(promoterDMR_df, deg, 
                          by.x = "nearestGene", by.y = "geneName"
                          )
  nonpromoterDMR_df <- merge(nonpromoterDMR_df, deg, 
                          by.x = "nearestGene", by.y = "geneName"
  )
  
  ### Summarise ###
  outFile <- sprintf("%s/dmrDelta_vs_degLog2FC_%s.pdf", outDir, dt)
  pdf(outFile, width = 19, height = 9)
  
  # subset by DEG FDR
  fdr_thresh = 0.05
  message(sprintf("Using %.2f as FDR threshold for DEGs", fdr_thresh))
  
  ## treat each DMR independently ##
  sub_promoterDMR_df <- promoterDMR_df[promoterDMR_df$FDR < fdr_thresh,]
  sub_nonpromoterDMR_df <- nonpromoterDMR_df[nonpromoterDMR_df$FDR < fdr_thresh,]
  
  plot_grid(plot_meth_rna_boxplot(sub_promoterDMR_df),
            plot_meth_rna_boxplot(sub_nonpromoterDMR_df),
            labels = c("Promoter DMR", "Non-promoter DMR"), label_size = 20)
  plot_meth_rna_scatter(promoterDMR_df, title = "Promoter DMR") + 
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  plot_meth_rna_scatter(sub_nonpromoterDMR_df, 
                        pearson_cor = F,
                        title = "Non-promoter DMR")
  
  dev.off()
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









