
library(gprofiler2)

#hypoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"
#hyperFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"

abcFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/Nasser-Neuronal-ABC_creTarget_hg38.bed"

rnaDEGfile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/231123/edgeR_RLVZvsSVZ_231123.txt"


hypoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/all_dmr_AME_fullMEME_251209/ame.tsv"

dmr2genes <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_link2Genes_ABC/251007/DMR_AnnotatedAll_251007.tsv"

#bgFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/MOTIFS_Hendrikse2022_RL_activeGenes.txt"
bgFile <- "/home/rstudio/isilon/src/transcription-factors/HOCOMOCOv12/H12CORE_meme_format.meme.motifs.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/TFBSenrichment"

outDir <- sprintf("%s/GOenrichment_%s", outDir, format(Sys.Date(),"%y%m%d"))
if (!file.exists(outDir)) dir.create(outDir, recursive = FALSE)

hypo <- read.delim(hypoFile, header = TRUE, sep = "\t", stringsAsFactors = 
FALSE)
##hyper <- read.delim(hyperFile, header = TRUE, sep = "\t", stringsAsFactors = 
##FALSE)  

# plot top 10 hits using light blue horizontal barplots. Plot the -log10 transofrm of the adj_p.value column. Use the motif_ID column for y-axis labels. Strip the suffix of evrything starting from the first dot.

hypo$motif_ID <- sub("\\..*", "", hypo$motif_ID)
##hyper$motif_ID <- sub("\\..*", "", hyper$motif_ID)  
hypo$log10adjp <- -log10(hypo$adj_p.value)
idx <- which(is.infinite(hypo$log10adjp))
hypo$log10adjp[idx] <-  max(hypo$log10adjp[!is.infinite(hypo$log10adjp)],na.rm=TRUE)+10
hypo <- hypo[!duplicated(hypo$motif_ID), ]  # Remove duplicates based on 
fg_genes <- hypo$motif_ID# [1:100]motif_ID
###hyper$log10adjp <- -log10(hyper$adj_p.value)
###idx <- which(is.infinite(hyper$log10adjp))
###hyper$log10adjp[idx] <- max(hyper$log10adjp[!is.infinite(hyper$log10adjp)],na.rm=TRUE)+10

##hyper <- hyper[!duplicated(hyper$motif_ID), ]  # Remove duplicates based on motif_IDsummary

bg <- read.delim(bgFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
bg_genes <- bg[,1]

dmr2genes_data <- read.delim(dmr2genes, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
fg_genes <- unique(dmr2genes_data$ABC_gene)

abc <- read.delim(abcFile, header=FALSE, stringsAsFactors = FALSE)
abc <-  subset(abc, V1 %in% paste("chr", c(1:22, "X", "Y"),sep=""))
bg_genes <- unique(abc[,4])

cat(sprintf("%i FG genes and %i BG genes loaded.\n", length(fg_genes), length(bg_genes)))

source("../../gprofiler2_helpers.R")
    # run pathway enrichment analysis (from `gprofiler2_helpers.R`)
    
    gost_res <- run_gost(
      query = fg_genes,
      organism = "hsapiens",
#      organism =  "gp__OUFl_gpIE_RyE",
      significant = TRUE,
      sources = c("GO:BP","GO:CC","REAC","WP","TF"),
      evcodes = TRUE,
      correction_method = "fdr",
      custom_bg = bg_genes,
      filename = file.path(outDir, "enriched_pathways")
    )

if (!is.null(gost_res) && nrow(gost_res$result) > 0) {

      res <- gost_res$result
      cat(sprintf("Found %i significantly enriched pathways.\n", nrow(res)))

      print(head(res$term_id))
      # save results table
      df <- gost_res$result
      df$parents <- NULL
      write.table(
        df,
        file = sprintf("%s/gost_enrichment_results.tsv", outDir),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    } else {
      cat("No significant enrichment found.\n")
    }

deg <- read.delim(rnaDEGfile, sep="\t", h=T, as.is=T)
