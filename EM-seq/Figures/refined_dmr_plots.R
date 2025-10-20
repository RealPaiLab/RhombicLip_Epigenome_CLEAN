library(dplyr)
library(ggplot2)
library(UpSetR)
library(glue)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../Analysis/FETHB2-FETHB3/utils.R")

### Env
dt <- format(Sys.Date(),"%y%m%d")
outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3",
               "/Figures/refined_dmr_plots/{dt}")

if (! dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}


### DMR histogram by DNAm difference ###
message("DMR histogram by DNAm difference")

## Configs
dmrFile <- get_configs("CPG_DMR_FILE")
outFile_png <- sprintf("%s/DMR_histogram_%s.png", outDir, dt)
outFile_pdf <- sprintf("%s/DMR_histogram_%s.pdf", outDir, dt)

## Load data
dmr_df <- read.table(dmrFile, stringsAsFactors = F, header = T)

## Pre-processing
dmr_df$status <- factor(ifelse(dmr_df$diff.Methy > 0, 'Hypo', 'Hyper'), 
                      levels = c('Hyper', 'Hypo'))
print(table(dmr_df$status))

## Plot
.plot <- ggplot(dmr_df, aes(x = diff.Methy*100, fill = status)) +
  geom_histogram(bins = round(sqrt(nrow(dmr_df))), color = "white") +
  theme_classic() +
  scale_fill_manual(values = c("Hypo" = "#756cb1", "Hyper" = "#141a86")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "DNAm Difference (%)", y = "Number of DMRs") + 
  annotate("segment", x = -5, xend = -60, y = 600, yend = 600, 
           color = "black", linewidth = 2,
           arrow = arrow(type = "closed", length = unit(0.2, "inches"))) +
  annotate("text", x = -30, y = 600, label = "in RL-VZ", vjust = -0.5, size = 8) +
  annotate("segment", x = 5, xend = 60, y = 600, yend = 600, 
           color = "black", linewidth = 2,
           arrow = arrow(type = "closed", length = unit(0.2, "inches"))) +
  annotate("text", x = 30, y = 600, label = "in RL-SVZ", vjust = -0.5, size = 8) +
  theme(legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20)
        )

## Saving
ggsave(outFile_png, plot = .plot, dpi = 600, height = 9, width = 10)
ggsave(outFile_pdf, plot = .plot, dpi = 600, height = 9, width = 10)



### DMR distribution by ENCODE elements ###
message("Generating DMR distribution by ENCODE elements")

## Configs
cpg_dmr_date <- get_configs("CPG_DMR_DATE")
cre_overlapFile <- glue("/.mounts/labs/pailab/private/projects",
                        "/FetalHindbrain/EMseq_FETHB3/output/downstream",
                        "/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded",
                        "/withoutBatchCorrection/{cpg_dmr_date}",
                        "/DMRs_cCRE_overlap_240712.txt")

gene_overlapFile <- glue("/.mounts/labs/pailab/private/projects",
                         "/FetalHindbrain/EMseq_FETHB3/output/downstream",
                         "/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/",
                         "withoutBatchCorrection/{cpg_dmr_date}",
                         "/DMRs_Gene_overlap_240712.txt")

outFile_png <- sprintf("%s/DMR_encodeIntersectUpset_%s.png", outDir, dt)
outFile_pdf <- sprintf("%s/DMR_encodeIntersectUpset_%s.pdf", outDir, dt)

## Load data
cre_overlap_df <- read.table(cre_overlapFile, stringsAsFactors = F, header = T)
gene_overlap_df <- read.table(gene_overlapFile, stringsAsFactors = F, header = T)

## Pre-processing
cre_overlap_df <- cre_overlap_df[, c("DMR.seqnames", "DMR.start", "DMR.end", "ENCODE.cCRE.Type")]
gene_overlap_df <- gene_overlap_df[! is.na(gene_overlap_df$GENE.seqnames),]
gene_overlap_df$ENCODE.cCRE.Type <- "geneBody"
gene_overlap_df <- gene_overlap_df[, c("DMR.seqnames", "DMR.start", "DMR.end", "ENCODE.cCRE.Type")]
overlap_df <- rbind(cre_overlap_df, gene_overlap_df)

creTypes <- overlap_df %>% 
  rename(cre = ENCODE.cCRE.Type) %>% 
  mutate(cre = stringr::str_remove(cre, ",CTCF-bound")) %>%
  select(cre) %>% 
  group_by(cre) %>% 
  summarise(unique(cre)) %>% 
  select(cre) %>% pull()

overlap_list <- overlap_df %>%
  mutate(dmr = paste(DMR.seqnames, DMR.start, DMR.end, sep = "_")) %>%
  rename(cre = ENCODE.cCRE.Type) %>%
  mutate(cre = stringr::str_remove(cre, ",CTCF-bound")) %>%
  select(dmr, cre) %>%
  group_by(cre) %>%
  group_split(.keep = T) %>%
  lapply(function(df) {df[[1]]}) %>%
  setNames(creTypes)

## Plot
png(outFile_png,width=7,height=5,units="in",res=600, pointsize=6)
upset(fromList(overlap_list), order.by = "freq", nsets = length(overlap_list))
dev.off()

pdf(outFile_pdf,width=7,height=5, pointsize = 6)
upset(fromList(overlap_list), order.by = "freq", nsets = length(overlap_list))
dev.off()



### DMR distribution by genomic annotations ###
message("Generating DMR distribution by genomic annotations")

## Configs
cpg_dmr_date <- get_configs("CPG_DMR_DATE")
feCB_overlapFile <- glue("/.mounts/labs/pailab/private/projects/",
                         "FetalHindbrain/EMseq_FETHB3/output/downstream/",
                         "EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded",
                         "/withoutBatchCorrection/{cpg_dmr_date}",
                         "/DMRs_feCB_overlap_240712.txt")

outFile_png <- sprintf("%s/DMR_fetalCerebellarIntersectUpset_%s.png", outDir, dt)
outFile_pdf <- sprintf("%s/DMR_fetalCerebellarIntersectUpset_%s.pdf", outDir, dt)

## Load data
feCB_overlap_df <- read.table(feCB_overlapFile, stringsAsFactors = F, header = T)

## Pre-processing
overlap_df <- feCB_overlap_df

creTypes <- overlap_df %>% 
  rename(cre = feCB_CRE.type) %>% 
  select(cre) %>% 
  group_by(cre) %>% 
  summarise(unique(cre)) %>% 
  select(cre) %>% pull()

overlap_list <- overlap_df %>%
  mutate(dmr = paste(DMR.seqnames, DMR.start, DMR.end, sep = "_")) %>%
  rename(cre = feCB_CRE.type) %>%
  select(dmr, cre) %>%
  group_by(cre) %>%
  group_split(.keep = T) %>%
  lapply(function(df) {df[[1]]}) %>%
  setNames(creTypes)

## Plot
png(outFile_png,width=7.5,height=5,units="in",res=600, pointsize=6)
upset(fromList(overlap_list), order.by = "freq", 
      nsets = length(overlap_list), 
      sets.bar.color = c("lightgrey", "#f8101b", "#f8101b", "#f8101b", "#fece2f"),
      text.scale = 1.75, point.size = 3,
      main.bar.color = c("lightgrey", rep("black", 13))
      )
dev.off()

pdf(outFile_pdf,width=7.5,height=5, pointsize=6)
upset(fromList(overlap_list), order.by = "freq", 
      nsets = length(overlap_list), 
      sets.bar.color = c("lightgrey", "#f8101b", "#f8101b", "#f8101b", "#fece2f"),
      text.scale = 1.75, point.size = 3,
      main.bar.color = c("lightgrey", rep("black", 13))
)
dev.off()






