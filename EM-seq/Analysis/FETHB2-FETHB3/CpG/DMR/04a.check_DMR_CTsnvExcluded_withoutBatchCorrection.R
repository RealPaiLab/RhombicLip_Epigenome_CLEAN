### NOTE: !!! remember to change date for DMRs.ENCODE.PLS.nearestTSS_240712.txt 
### and DMRs.ENCODE.ELS.nearestTSS_240712.txt if you rerun anything. These cannot
### be fixed because they depends on when you run the previous annotation step

rm(list=ls())

library(glue)
library(methylKit)
library(dplyr)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")

dt <- format(Sys.Date(),"%y%m%d")
cpg_dmr_date <- get_configs("CPG_DMR_DATE")

rootDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")

outDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection/%s/checkDMR",
                  rootDir, 
                  cpg_dmr_date
                  )

logFile <- sprintf("%s/check_DMR_CTsnvExcluded_withoutBatchCorrection_%s.log",
                   outDir, 
                   dt
)

main <- function() {
  ### Load data ###
  ## Hendrikse RL-VZ/SVZ top 50 DEGs (from glutamatergic lineage snRNA-seq)
  vz_deg <-readxl::read_excel(
    glue("/.mounts/labs/pailab/src/gene_list/brain-development",
         "/Hendrikse2022_snRNA_glutamatergicTopDEGs.xlsx"
    ),
    sheet = 2
  )
  svz_deg <-readxl::read_excel(
    glue("/.mounts/labs/pailab/src/gene_list/brain-development",
         "/Hendrikse2022_snRNA_glutamatergicTopDEGs.xlsx"
    ),
    sheet = 3
  )
  
  ## Annotated DMRs ##
  dmr_ol_els <- read.table(
    glue("/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3",
         "/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
         "/CTsnv_excluded/withoutBatchCorrection/{cpg_dmr_date}",
         "/DMRs.ENCODE.ELS.nearestTSS_240712.txt"
    ),
    stringsAsFactors = F, 
    header = T
  )
  
  dmr_ol_pls <- read.table(
    glue("/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3",
         "/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
         "/CTsnv_excluded/withoutBatchCorrection/{cpg_dmr_date}",
         "/DMRs.ENCODE.PLS.nearestTSS_240712.txt"
    ),
    stringsAsFactors = F, 
    header = T
  )
  
  ## methylation data ##
  # load regport #
  report_dir <- get_configs("CPG_REPORT_DIR")
  
  files <- list.files(
    report_dir,
    pattern = "*.txt.gz$")
  
  sample.id <- extract_sampleId(files)
  
  message("Loading cytosine reports")
  myobj <- methRead(as.list(paste(report_dir, files, sep = "/")),
                    sample.id = as.list(sample.id),
                    assembly = "hg38",
                    treatment = get_TissueType(sample.id),
                    context = "CpG", 
                    pipeline = "bismarkCytosineReport",
                    mincov = 1 # use mincov 1 as CpGs will be grouped by windows
  )
  
  # subset by interested regions #
  interested_regions <- rbind(dmr_ol_els, dmr_ol_pls) %>%
    dplyr::select(c(1,2,3)) %>% 
    dplyr::rename(seqnames = DMR.seqnames, start = DMR.start, end = DMR.end) %>%
    distinct()
  
  interested_gr <- GenomicRanges::makeGRangesFromDataFrame(interested_regions)
  subobj <- selectByOverlap(myobj, interested_gr)
  rm(myobj); gc()
  
  
  
  ### Utils ###
  
  #' Summarize sample percent methylated Cs of the given region
  #' @param meth (methylRawList) A methylRawList object generated using MethylKit
  #' @param chr (character) Chromosome of the interested region (e.g. "chr1")
  #' @param start (numeric) Start position of the interested region
  #' @param end (numeric) End position of the interested region
  #' @param return (data.frame)
  summarize_region_meth <- function(meth, chr, start, end) {
    tmp <- selectByOverlap(meth, 
                           GRanges(seqnames = Rle(chr),
                                   ranges = IRanges(start = as.integer(start), 
                                                    end = as.integer(end)
                                   )
                           )
    )
    united_meth <- getData(unite(tmp, destrand=T, min.per.group = 1L))
    
    c_cols <- which(grepl("numCs", colnames(united_meth)))
    cov_cols <- which(grepl("coverage", colnames(united_meth)))
    
    res_df <- data.frame(
      Meth_c = colSums(united_meth[,c_cols]),
      Coverage = colSums(united_meth[,cov_cols]),
      Percent = colSums(united_meth[,c_cols])/colSums(united_meth[,cov_cols])*100,
      Tissue = factor(ifelse(meth@treatment, "VZ", "SVZ"), levels = c("VZ", "SVZ"))
    )
    
    return(res_df)
  }
  
  
  #' Plot sample percent methylated Cs of the given region grouped by tissue type
  #' @param meth (methylRawList) A methylRawList object generated using MethylKit
  #' @param chr (character) Chromosome of the interested region (e.g. "chr1")
  #' @param start (numeric) Start position of the interested region
  #' @param end (numeric) End position of the interested region
  #' @param title (character) The title of the output plot
  #' @param return (ggplot) The violin ggplot object 
  plot_region_meth <- function(meth, chr, start, end, title) {
    message(glue("Start to process for: {title}"))
    
    res_df <- summarize_region_meth(meth, chr, start, end)
    
    .plot <- ggplot(res_df, aes(x = Tissue, y = Percent, fill = Tissue)) + 
      geom_violin(color = "NA") + 
      scale_fill_manual(values = c("VZ" = "gold1", "SVZ" = "forestgreen")) +
      geom_jitter(aes(color = Coverage)) +
      scale_color_gradient(low = "lightgrey", high = "black") +
      ggtitle(title) + 
      ylab("Percent Methylated Cs (%)") +
      theme_classic() + 
      theme(axis.text = element_text(size = 12), 
            axis.title = element_text(size = 12)
      ) +
      ylim(0,100)
    
    return(.plot)
  }
  
  
  #' Plot sample percent methylated Cs of the dmr_ol dataframe grouped by tissue type 
  #' @param additionalTitle (character) Addtional string goes into the plot title
  #' @return (list) A list of ggplot objects; one for each gene.
  plot_dmr_ol_df <- function(df, additionalTitle = "") {
    plot_list <- apply(df, 
                       1,  
                       function(x) {
                         plot_region_meth(subobj, 
                                          x["DMR.seqnames"], 
                                          x["DMR.start"], 
                                          x["DMR.end"], 
                                          sprintf("%s%s (areaStat = %.2f)",
                                                  additionalTitle,
                                                  x["nearestTSS_proteinCoding"], 
                                                  as.numeric(x["DMR.areaStat"])
                                          )
                         )
                       }
    )
    
    #recordPlot(gridExtra::grid.arrange(grobs = plot_list))
    .plot <- gridExtra::arrangeGrob(grobs = plot_list)
    #.plot <- recordPlot()
    return(.plot)
  }
  
  
  ### Check DMRs near RL-VZ/SVZ DEGs ###
  ## pls ##
  # VZ #
  dmr_ol_pls_vz <- merge(dmr_ol_pls, 
                         vz_deg, 
                         by.x = "nearestTSS_proteinCoding",
                         by.y = "Gene"
  )
  message(sprintf("%i PLS DMRs near RL-VZ DEGs", nrow(dmr_ol_pls_vz)))
  message(sprintf("%i DMRs with areaStat < 0, %i DMR with areaStat > 0",
                  sum(dmr_ol_pls_vz$DMR.areaStat < 0),
                  sum(dmr_ol_pls_vz$DMR.areaStat > 0)
  )
  )
  
  plot_dmr_ol_df(dmr_ol_pls_vz,
                 additionalTitle = "DMR ol PLS near RL-VZ gene\n")
  outFile <- glue("{outDir}/dmr_ol_pls_vz_violins_{dt}.png")
  ggsave(outFile, dpi = 600, width = 7, height = 7)
  outFile <- glue("{outDir}/dmr_ol_pls_vz_violins_{dt}.pdf")
  ggsave(outFile, dpi = 600, width = 7, height = 7)
  
  # SVZ #
  dmr_ol_pls_svz <- merge(dmr_ol_pls, 
                          svz_deg, 
                          by.x = "nearestTSS_proteinCoding",
                          by.y = "Gene"
  )
  message(sprintf("%i PLS DMRs near RL-SVZ DEGs", nrow(dmr_ol_pls_svz)))
  message(sprintf("%i DMRs with areaStat < 0, %i DMR with areaStat > 0",
                  sum(dmr_ol_pls_svz$DMR.areaStat < 0),
                  sum(dmr_ol_pls_svz$DMR.areaStat > 0)
  )
  )
  
  .plot <- plot_dmr_ol_df(dmr_ol_pls_svz,
                 additionalTitle = "DMR ol PLS near RL-SVZ gene\n")
  outFile <- glue("{outDir}/dmr_ol_pls_svz_violins_{dt}.png")
  ggsave(outFile, plot = .plot, dpi = 600, width = 10, height = 10)
  outFile <- glue("{outDir}/dmr_ol_pls_svz_violins_{dt}.pdf")
  ggsave(outFile, plot = .plot, dpi = 600, width = 10, height = 10)
  
  ## els ##
  # VZ #
  dmr_ol_els_vz <- dmr_ol_els[dmr_ol_els$nearestTSS_proteinCoding %in% vz_deg$Gene,]
  
  dmr_ol_els_vz <- merge(dmr_ol_els, 
                         vz_deg, 
                         by.x = "nearestTSS_proteinCoding",
                         by.y = "Gene"
  )
  message(sprintf("%i ELS DMRs near RL-VZ DEGs", nrow(dmr_ol_els_vz)))
  message(sprintf("%i DMRs with areaStat < 0, %i DMR with areaStat > 0",
                  sum(dmr_ol_els_vz$DMR.areaStat < 0),
                  sum(dmr_ol_els_vz$DMR.areaStat > 0)
                  )
          )
  
  # SVZ #
  dmr_ol_els_svz <- merge(dmr_ol_els, 
                          svz_deg, 
                          by.x = "nearestTSS_proteinCoding",
                          by.y = "Gene"
                          )
  message(sprintf("%i ELS DMRs near RL-SVZ DEGs", nrow(dmr_ol_els_svz)))
  message(sprintf("%i DMRs with areaStat < 0, %i DMR with areaStat > 0",
                  sum(dmr_ol_els_svz$DMR.areaStat < 0),
                  sum(dmr_ol_els_svz$DMR.areaStat > 0)
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









