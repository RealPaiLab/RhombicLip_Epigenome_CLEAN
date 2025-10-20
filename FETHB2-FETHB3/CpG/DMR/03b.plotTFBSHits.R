# Plot top hits for TFBS enrichment analysis
# for differentially methylated regions (DMRs) in FETHB2-FETHB3 dataset
library(ggplot2)

hypoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"
hyperFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/TFBSenrichment"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}


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

hypoTop <- hypo[1:10, ]
#hypoTop$motif_ID <- factor(hypoTop$motif_ID, levels = hypoTop$motif_ID[order(hypoTop$log10adjp, decreasing = FALSE)])
hyperTop <- hyper[1:10, ]
#hyperTop$motif_ID <- factor(hyperTop$motif_ID, levels = hyperTop$motif_ID[order(hyperTop$log10adjp, decreasing = FALSE)])

ggplot(hypoTop, aes(x = reorder(motif_ID, log10adjp), y = log10adjp)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  labs(title = "Top 10 Hypomethylated DMRs TFBS Hits",
       x = "Motif ID",
       y = "-log10(adjusted p)") +
  theme_minimal(base_size = 24)
ggsave(sprintf("%s/hypo_top_hits.pdf", outDir), width = 8, height = 6)

ggplot(hyperTop, aes(x = reorder(motif_ID, log10adjp), y = log10adjp)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  labs(title = "Top 10 Hypermethylated DMRs TFBS Hits",
       x = "Motif ID",
       y = "-log10(adjusted p)") +
  theme_minimal(base_size = 24) 
ggsave(sprintf("%s/hyper_top_hits.pdf", outDir), width = 8, height = 6)  

# combine hyper and hypo into one plot. Colour the hypo bars with #7B3294 and hyper bars with #008837
combined <- rbind(hypoTop, hyperTop)
combined$Type <- rep(c("Hypomethylated", "Hypermethylated"), each = 10)

ggplot(combined, aes(x = reorder(motif_ID, log10adjp), y = log10adjp, fill = Type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Hypomethylated" = "#7B3294", "Hypermethylated" = "#008837")) +
  labs(title = "Top 10 DMRs TFBS Hits",
       x = "Motif ID",
       y = "-log10(adjusted p)") +
  theme_minimal(base_size = 24)
ggsave(sprintf("%s/combined_top_hits.pdf", outDir), width = 8, height = 6)
