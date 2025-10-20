rm(list=ls())

library(GenomicRanges)
library(dplyr)
library(glue)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")

### dir config ###
rootDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")

segDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG",
               "/DMVs/CTsnv_excluded",
               "/240712",
               "/Segmentation"
               )
dmrFile <- get_configs("CPG_DMR_FILE")
geneFile <- get_configs("GENCODE_GENE_FILE")

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/DMVs/CTsnv_excluded/%s",rootDir, "240712")
logFile <- sprintf("%s/DNAmValley_%s.log",outDir, dt)


main <- function() {
  message(sprintf("Output will be written under: %s\n",outDir))

  ### Load DMR ###
  ranges_dmrs <- get_cpg_dmrs()
  
  ### Load segmentation ###
  message(sprintf("Loading segmentation files from %s", segDir))
  segFiles <- list.files(segDir, pattern = "FETHB[0-9]_.+_seg.bed")
  
  mv_list <- NULL
  for (segFile in segFiles) {
    sample <- substr(segFile, 1, 19)
    message(sprintf("Processing %s", segFile))
    
    seg <- read.table(sprintf("%s/%s", segDir, segFile), skip = 1)
    colnames(seg) <- c("chr", "start", "end", 
                       "segCluster", "meanCpGm", "strand", 
                       "V7", "V8", 
                       "RGB"
                       ) #cluster1 is lowest
    
    message(sprintf("Subsetting DNAm Valleys for: %s", sample))
    # mv definition, 0.15 CpG methylation and over 5kb length; 
    # https://doi.org/10.1186/s13059-018-1390-8 
    mv <- seg %>%
      mutate(width = end - start) %>%
      filter(width >= 5000 & meanCpGm <= 15) %>% 
      dplyr::select(c("chr", "start", "end", "segCluster", "meanCpGm", "width"))
    message(
      sprintf(
        "%d DNAm Valleys identified, with median length %.2f and median meanCpG methyaltion %.2f\n", 
        nrow(mv), 
        median(mv$width), 
        median(mv$meanCpGm)
        )
      )
    
    mv_list[[sample]] <- mv
  }
  
  
  ### MV Gene overlap ###
  message("Overlapping MVs with GENCODE gene annotation")
  
  ## Load GENCODE db ##
  genes <- rtracklayer::readGFF(geneFile)
  
  # GENCODE GRange Object - all gene_type #
  geneGR_all <- get_gencode_anno(gene_type = "all")
  
  # GENCODE GRange Object - only protein coding #
  geneGR_proteinCoding <- get_gencode_anno(gene_type = "protein_coding")
  
  ## annotate ##
  GR_use <- geneGR_proteinCoding
  #GR_use <- geneGR_all
  for (i in 1:length(mv_list)) {
    sample <- names(mv_list)[[i]]
    mv <- mv_list[[i]]
    message(sprintf("Processing %s", sample))
    
    ranges_mv <- GRanges(mv$chr, 
                         IRanges(mv$start, mv$end))
    mcols(ranges_mv)$segCluster <- mv$segCluster
    mcols(ranges_mv)$meanCpGm <- mv$meanCpGm
    
    
    ol <- findOverlaps(
      ranges_mv, GR_use, minoverlap = 0
    )
    uq_query <- length(unique(queryHits(ol)))
    uq_sbj <- length(unique(subjectHits(ol)))
    
    message(sprintf("Found total %i OL; %i MV X %i Genes", 
                    length(ol), uq_query, uq_sbj
    ))
    
    # construct df #
    qry <- as.data.frame(ranges_mv[queryHits(ol)])
    qry_noHit <- as.data.frame(ranges_mv[-queryHits(ol)])
    
    sbj <- as.data.frame(GR_use[subjectHits(ol)])
    colnames(sbj) <- paste("GENE",colnames(sbj),sep=".")
    
    regions <- cbind(qry,sbj)
    
    tmp <- cbind(qry_noHit, 
                 data.frame(matrix(nrow = nrow(qry_noHit), ncol = ncol(sbj)))
                 )
    colnames(tmp) <- colnames(regions)
    tmp$GENE.name <- "noMatch"
    
    regions <- rbind(regions, tmp)
    
    # simplify #
    simp_regions <- regions %>% 
      group_by(seqnames, start, end, width, segCluster, meanCpGm) %>% 
      summarise(
        GENE.name = toString(sort(unique(GENE.name))), 
        n_gene_overlap = ifelse(GENE.name == "noMatch", 0, n()), 
      ) %>%
      arrange(.by_group = T)
    
    # update mv_list #
    colnames(simp_regions)[which(colnames(simp_regions) == "seqnames") ] <- "chr"
    mv_list[[sample]] <- simp_regions 
  }
  
  
  ### Mark MVs with DMR overlaps ###
  message("Marking MVs with DMR overlaps")
  dmv_list <- NULL
  for (i in 1:length(mv_list)) {
    sample <- names(mv_list)[[i]]
    mv <- mv_list[[i]]
    message(sprintf("Processing %s", sample))
    
    ranges_mv <- GRanges(mv$chr, 
                           IRanges(mv$start, mv$end))
    mcols(ranges_mv)$segCluster <- mv$segCluster
    mcols(ranges_mv)$meanCpGm <- mv$meanCpGm
    mcols(ranges_mv)$GENE.name <- mv$GENE.name
    mcols(ranges_mv)$n_gene_overlap <- mv$n_gene_overlap
  
    # findOverlaps can't work with percentage overlap; 
    # will need to switch to bedtools, using 25 which is half of smallest DMR
    ol <- findOverlaps(
      ranges_dmrs, ranges_mv, minoverlap = 25 
    )
    uq_query <- length(unique(queryHits(ol)))
    uq_sbj <- length(unique(subjectHits(ol)))
    message(sprintf("Found total %i OL; %i DMR X %i MV", 
                    length(ol), uq_query, uq_sbj
    ))
    
    qry <- as.data.frame(ranges_dmrs[queryHits(ol)])
    colnames(qry) <- paste("DMR",colnames(qry),sep=".")
    qry_noHit <- as.data.frame(ranges_dmrs[-queryHits(ol)])
    
    sbj <- as.data.frame(ranges_mv[subjectHits(ol)])
  
    dmv <- cbind(qry,sbj) %>%
      mutate(svzDNAm = ifelse(DMR.diff.Methy > 0, "Hypo", "Hyper")) %>%
      group_by(seqnames, start, end,
               width, segCluster, meanCpGm, 
               GENE.name, n_gene_overlap
               ) %>%
      summarise(n_DMR_overlap = n(), 
                DMRs = toString(
                  sort(sprintf("%s:%s-%s", DMR.seqnames, DMR.start, DMR.end))),
                svzDNAm_hypo_dmr = sum(svzDNAm == "Hypo"), 
                svzDNAm_hyper_dmr = sum(svzDNAm == "Hyper"),
                svzDNAm = toString(sort(unique(svzDNAm)))
                ) #%>% filter(svzDNAm != "Hyper, Hypo") # remove ambiguous DMVs
    
    message(
      sprintf("After filtering, %d dmv marked, %d are hypomethylated in RL-SVZ, %d are hypermethylated.\n", 
              nrow(dmv), 
              sum(dmv$svzDNAm == "Hypo"), 
              sum(dmv$svzDNAm == "Hyper")
              )
            )
    
    colnames(dmv)[which(colnames(dmv) == "seqnames") ] <- "chr"
    dmv$region <- ifelse(get_TissueType(sample) == 1, "VZ", "SVZ")
  
    dmv$sample <- sample
    dmv_list[[sample]] <- dmv
  }
  
  ## Count DMV overlapping genes across samples
  c <- sort(
    table(
      unlist(
        lapply(dmv_list, 
               function(x) {paste(x$region, x$svzDNAm, x$GENE.name, sep = "_")}
               )
        )
      )
    )/length(dmv_list)
  
  c <- c[! grepl("noMatch", names(c))]
  c <- as.data.frame(c)
  c$region <- stringr::str_split(c$Var1, "_", simplify = T)[,1]
  c$meth <- stringr::str_split(c$Var1, "_", simplify = T)[,2]
  c$gene <- stringr::str_split(c$Var1, "_", simplify = T)[,3]
  orders <- c %>%
    group_by(gene) %>%
    summarise(freq = sum(Freq)) %>%
    arrange(freq) %>%
    pull(gene)
  c$gene <- factor(c$gene, levels = orders)
  
  
  ggplot(as.data.frame(c), 
         aes(x = gene, y = ifelse(region == "SVZ", Freq, -Freq), fill = meth)
         ) + 
    geom_bar(stat = "identity", position = "identity") + 
    scale_y_continuous(labels = abs) +
    coord_flip() + 
    theme_bw() +
    ylab("Percent sample") +
    xlab("Overlapping gene(s)") + 
    ggtitle("VZ vs SVZ") +
    theme(axis.text.y = element_text(size = 7)) +
    geom_hline(yintercept = 0)
    
  ### Merge sample DMVs ###
  all_dmv <- do.call(rbind, dmv_list)
  ranges_dmv <- makeGRangesFromDataFrame(all_dmv, keep.extra.columns = TRUE)
  ranges_dmv_merged <- reduce(ranges_dmv)
  
  ol <- findOverlaps(
    ranges_dmv, ranges_dmv_merged, minoverlap = 0
  )
  
  merged_range_df <- as.data.frame(ranges_dmv_merged)[subjectHits(ol),]
  merged_range_df <- merged_range_df[,c("start", "end", "width")]
  colnames(merged_range_df) <- paste("merged", 
                                     colnames(merged_range_df), 
                                     sep = "_"
                                     )
  all_dmv <- cbind(all_dmv, merged_range_df)
  
  message(sprintf("%d unique DMVs (median length %d; %d - %d) found after merging all %d DMVs (median length %d; %d - %d)", 
                  length(ranges_dmv_merged), 
                  median(as.data.frame(ranges_dmv_merged)$width), 
                  min(all_dmv$merged_width),
                  max(all_dmv$merged_width),
                  nrow(all_dmv), 
                  median(all_dmv$width),
                  min(all_dmv$width),
                  max(all_dmv$width)
                  )
          )
  
  outFile <- sprintf("%s/DMVs_%s.tsv", outDir, dt)
  write.table(all_dmv, outFile, col.names = T, row.names = F, quote = F, sep = "\t")
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


















