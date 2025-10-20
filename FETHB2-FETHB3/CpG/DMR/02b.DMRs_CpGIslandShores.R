# annotate DMRs by CpG island, shores, etc.,

library(annotatr)
library(GenomicRanges)
# make a pie chart with the annot.type
library(ggplot2)
library(dplyr)

CPG_DMR_FILE <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_CpGannot"
dt <- format(Sys.time(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
}

logFile <- sprintf("%s/DMRs_CpGIslandShores.log", outDir)
sink(logFile,split=TRUE)
tryCatch({
  
dmr <- read.delim(CPG_DMR_FILE, stringsAsFactors = FALSE,sep="\t",header=TRUE)

dmrGR <- GRanges(
  seqnames = dmr$chr,
  ranges = IRanges(start = dmr$start, end = dmr$end)
)

annots <- c("hg38_cpg_islands",
            "hg38_cpg_shores",
            "hg38_cpg_shelves",
            "hg38_cpg_inter")

annotations <- build_annotations(annotations = annots,
                                  genome = "hg38")

dm_annotated <- annotate_regions(
  regions = dmrGR,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE
)                                

df_dm_annotated <- as.data.frame(dm_annotated)

cat("Distribution of DMRs by CpG annotation:\n")
print(table(df_dm_annotated$annot.type))


df_dm_annotated$annot.type <- factor(df_dm_annotated$annot.type,
                                     levels = c("hg38_cpg_islands",
                                                "hg38_cpg_shores",
                                                "hg38_cpg_shelves",
                                                "hg38_cpg_inter"))

p <- df_dm_annotated %>%
  group_by(annot.type) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = "", y = count, fill = annot.type)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c(
    "hg38_cpg_islands" = "darkolivegreen3",
    "hg38_cpg_shores" = "chocolate4",
    "hg38_cpg_shelves" = "cornflowerblue",
    "hg38_cpg_inter" = "blue4"
  )) +  
    coord_polar(theta = "y") +
  theme_void()

# add labels to the pie chart with the percentage of each type
p <- p + geom_text(aes(label = scales::percent(count / sum(count))),
                   position = position_stack(vjust = 0.5), color = "white",
                   size = 10) 
ggsave(p,
    file=sprintf("%s/DMRs_CpGIslandShores_pie.pdf", outDir),
    width = 8, height = 8)

}, error=function(ex){
    print(ex)
}, finally={
    sink(NULL)
    print("Done")
  }
)
