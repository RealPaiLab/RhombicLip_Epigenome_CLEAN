
library(gprofiler2)

hypoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"
hyperFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"

bgFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/MOTIFS_Hendrikse2022_RL_activeGenes.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/TFBSenrichment"

outDir <- sprintf("%s/GOenrichment_%s", outDir, format(Sys.Date(),"%y%m%d"))
if (!file.exists(outDir)) dir.create(outDir, recursive = FALSE)

hypo <- read.delim(hypoFile, header = TRUE, sep = "\t", stringsAsFactors = 
FALSE)
hyper <- read.delim(hyperFile, header = TRUE, sep = "\t", stringsAsFactors = 
FALSE)  

# plot top 10 hits using light blue horizontal barplots. Plot the -log10 transofrm of the adj_p.value column. Use the motif_ID column for y-axis labels. Strip the suffix of evrything starting from the first dot.

# replace motif name NDF1 with NEUROD1
hypo$motif_ID <- gsub("NDF1", "NEUROD1", hypo$motif_ID)
hyper$motif_ID <- gsub("NDF1", "NEUROD1", hyper$motif_ID)
# replace motif name NDF2 with NEUROD2
hypo$motif_ID <- gsub("NDF2", "NEUROD2", hypo$motif_ID)
hyper$motif_ID <- gsub("NDF2", "NEUROD2", hyper$motif_ID)

hypo$motif_ID <- sub("\\..*", "", hypo$motif_ID)
hyper$motif_ID <- sub("\\..*", "", hyper$motif_ID)  
hypo$log10adjp <- -log10(hypo$adj_p.value)
idx <- which(is.infinite(hypo$log10adjp))
hypo$log10adjp[idx] <-  max(hypo$log10adjp[!is.infinite(hypo$log10adjp)],na.rm=TRUE)+10
hypo <- hypo[!duplicated(hypo$motif_ID), ]  # Remove duplicates based on motif_ID
hyper$log10adjp <- -log10(hyper$adj_p.value)
idx <- which(is.infinite(hyper$log10adjp))
hyper$log10adjp[idx] <- max(hyper$log10adjp[!is.infinite(hyper$log10adjp)],na.rm=TRUE)+10

hyper <- hyper[!duplicated(hyper$motif_ID), ]  # Remove duplicates based on motif_ID



bg <- read.delim(bgFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
bg[,1] <- sub("MOTIF ", "",bg[,1])
dpos <- regexpr("\\.", bg[,1])
bg_genes <- substr(bg[,1], 1, dpos-1)

source("../../gprofiler2_helpers.R")
    # run pathway enrichment analysis (from `gprofiler2_helpers.R`)
    fg_genes <- hypo$motif_ID
    gost_res <- run_gost(
      query = fg_genes,
      organism = "gp__OUFl_gpIE_RyE",
      significant = FALSE,
      evcodes = TRUE,
      correction_method = "fdr",
      custom_bg = bg_genes,
      filename = file.path(outDir, "enriched_pathways")
    )
