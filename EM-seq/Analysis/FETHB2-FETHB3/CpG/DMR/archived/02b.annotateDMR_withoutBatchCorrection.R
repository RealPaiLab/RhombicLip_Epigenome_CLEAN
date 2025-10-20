# adapted from https://github.com/RealPaiLab/FetalHindbrain_Epigenetics/blob/master/FET_HB2/dmrAnnotate.R

## ---- Identification of ccREs which overlap with DMRs

rm(list=ls())

library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(AnnotationHub)
library(UCSCRepeatMasker)
library(reshape2)
library(networkD3)


### config ###
dt <- format(Sys.Date(),"%y%m%d")

rootDir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG"
outDir <- sprintf("%s/DMRs/CTsnv_included/withoutBatchCorrection/%s",rootDir, "240419")

ccre_dir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/GRCh38-cCREs.bed"
dmrDir <- outDir

geneFile <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/gencode.v42.basic.annotation.gtf"

teFile <- "/.mounts/labs/pailab/private/xsun/Database/RepeatMasker/RepeatMasker_open406_Dec2013_Dfam20_hg38/hg38.fa.out.gz"

aldingerCB_CRE_dir <- "/.mounts/labs/pailab/private/xsun/output/ncMutMB/20240314/cre/Aldinger-FetalCB"

logFile <- sprintf("%s/annotateDMR_dssBatchCorrection_%s.log",outDir, dt)

main <- function() {
  message(sprintf("Output will be written under: %s\n",outDir))
  message(sprintf("cCRE input: %s", ccre_dir))
  message(sprintf("DMR input: %s", dmrDir))
  
  ### read files ###
  ccres <- read.delim(
    ccre_dir,
    header=F,sep="\t"
  )
  
  message("Intersecting DMRs with cCREs...")
  
  dmrs <- read.delim(
    sprintf("%s/DMRs.csv",dmrDir),
    sep="\t")
  
  
  ### cCRE overlap ###
  df_dmrs <- dmrs
  df_ccres <- ccres
  
  ranges_ccres <- GRanges(ccres$V1, IRanges(ccres$V2, ccres$V3))
  ranges_ccres$ID <- ccres$V5
  ranges_ccres$Type <- ccres$V6
  
  ranges_dmrs <- GRanges(df_dmrs$chr, 
                         IRanges(df_dmrs$start, df_dmrs$end))
  mcols(ranges_dmrs)$nCG <- df_dmrs$nCG
  mcols(ranges_dmrs)$areaStat <- df_dmrs$areaStat
  
  ## Check overlap
  ol <- findOverlaps(
    ranges_dmrs, ranges_ccres, minoverlap = 0
  )
  
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  
  message(sprintf("Found total %i OL; %i DMR X %i cCRE", 
                  length(ol), uq_query, uq_sbj
  ))
  
  ## Summarise cCRE overlap 
  qry <- as.data.frame(ranges_dmrs[queryHits(ol)])
  colnames(qry) <- paste("DMR",colnames(qry),sep=".")
  qry_noHit <- as.data.frame(ranges_dmrs[-queryHits(ol)])
  
  sbj <- as.data.frame(ranges_ccres[subjectHits(ol)])
  colnames(sbj) <- paste("ENCODE.cCRE",colnames(sbj),sep=".")
  
  regions <- cbind(qry,sbj)
  
  tmp <- cbind(qry_noHit, data.frame(matrix(nrow = nrow(qry_noHit), ncol = ncol(sbj))))
  colnames(tmp) <- colnames(regions)
  tmp$ENCODE.cCRE.Type <- "noMatch"
  
  regions <- rbind(regions, tmp)
  
  simp_regions <- regions %>% 
    mutate(simplifiedType = sub(",CTCF-bound", "",ENCODE.cCRE.Type)) %>% # simplify by removing CTCF-bound term
    group_by(DMR.seqnames, DMR.start, DMR.end, DMR.width, DMR.strand, DMR.nCG, DMR.areaStat) %>% 
    summarise(
      ENCODE.cCRE.Type = toString(sort(unique(simplifiedType))), 
      ENCODE.cCRE.N_overlap = ifelse(ENCODE.cCRE.Type == "noMatch", 0, n())
    ) %>% 
    arrange(.by_group = T)
  
  ## save cCRE overlap regions and simplified regions 
  write.table(regions,
              file=sprintf("%s/DMRs_cCRE_overlap_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  
  write.table(simp_regions,
              file=sprintf("%s/DMRs_cCRE_overlap_simplified_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  rm(ol,qry,sbj)
  
  ### transposable elements
  message("Intersecting DMRs with TEs...")
  
  ## Load TE db
  ah <- AnnotationHub()
  query(ah, c("RepeatMasker", "Homo sapiens"))
  rmskhg38 <- ah[["AH99003"]]
  
  ## Check overlap
  ol <- findOverlaps(
    ranges_dmrs, rmskhg38, minoverlap = 0
  )
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  
  message(sprintf("Found total %i OL; %i DMR X %i TEs", 
                  length(ol), uq_query, uq_sbj
  ))
  
  
  ## summarise TE count
  qry <- as.data.frame(ranges_dmrs[queryHits(ol)])
  colnames(qry) <- paste("DMR",colnames(qry),sep=".")
  qry_noHit <- as.data.frame(ranges_dmrs[-queryHits(ol)])
  
  sbj <- as.data.frame(rmskhg38[subjectHits(ol)])
  colnames(sbj) <- paste("TE",colnames(sbj),sep=".")
  
  TEregions <- cbind(qry,sbj)
  
  tmp <- cbind(qry_noHit, data.frame(matrix(nrow = nrow(qry_noHit), ncol = ncol(sbj))))
  colnames(tmp) <- colnames(TEregions)
  tmp$TE.repName <- "noMatch"
  tmp$TE.repClass <- "noMatch"
  tmp$TE.repFamily <- "noMatch"
  
  TEregions <- rbind(TEregions, tmp)
  
  simp_TEregions <- TEregions %>% 
    group_by(DMR.seqnames, DMR.start, DMR.end, DMR.width, DMR.strand, DMR.nCG, DMR.areaStat) %>% 
    summarise(
      TE.repName = toString(sort(unique(TE.repName))), 
      TE.repClass = toString(sort(unique(TE.repClass))), 
      TE.repFamily = toString(sort(unique(TE.repFamily))), 
      TE.N_overlap = ifelse(TE.repClass == "noMatch", 0, n()), 
      TE.ambiguous = all(grepl("\\?", TE.repClass))
    ) %>%
    arrange(.by_group = T)
  
  ## Save TE overlap regions and simplified regions
  write.table(TEregions,
              file=sprintf("%s/DMRs_TE_overlap_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  
  write.table(simp_TEregions,
              file=sprintf("%s/DMRs_TE_overlap_simplified_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  rm(ol,qry,sbj)
  
  
  ### Count Enhancer/Promoter Overlap
  ## for enhancers
  message("Enhancer overlap ***** ")
  idx <- which(ranges_ccres$Type %in% 
                 c("dELS,CTCF-bound","pELS","pELS,CTCF-bound","dELS"))
  enh <- ranges_ccres[idx]
  message(sprintf("%i enhancers",length(enh)))
  enh  <- reduce(enh) # bedtools merge
  message(sprintf("%i reduced ranges",length(enh)))
  
  message(sprintf("Input: %s DMRs, %s ENCODE cCREs",
                  prettyNum(length(ranges_dmrs),big.mark=","), 
                  prettyNum(length(ranges_ccres),big.mark=",")
  ))
  ol <- findOverlaps(
    ranges_dmrs, enh
  )
  
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  
  message(sprintf("Found total %i OL; %i DMR (%1.2f %%) overlap %i enhancers", 
                  length(ol), uq_query, 
                  round((uq_query/length(ranges_dmrs))*100,1),
                  uq_sbj
  ))
  rm(enh)
  
  ## for promoters
  message("Promoter overlap ***** ")
  idx <- which(ranges_ccres$Type %in% 
                 c("PLS,CTCF-bound","PLS"))
  prom <- ranges_ccres[idx]
  message(sprintf("%i promoters",length(prom)))
  prom  <- reduce(prom)
  message(sprintf("%i reduced ranges",length(prom)))
  
  ol <- findOverlaps(
    ranges_dmrs, prom
  )
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  message(sprintf("Found total %i OL; %i DMR (%1.2f %%) overlap %i promoters", 
                  length(ol), uq_query, 
                  round((uq_query/length(ranges_dmrs))*100,1),
                  uq_sbj
  ))
  rm(ol,uq_query,uq_sbj)
  
  
  ### generate summary plot
  ## complete cCRE
  x <- table(regions$ENCODE.cCRE.Type)
  y <- data.frame(
    ENCODE.cCRE.Category=names(x),
    numDMRs_overlapping=as.numeric(x),
    pctOverlap=(as.numeric(x)/nrow(regions))
  )
  y <- y[order(y$numDMRs_overlapping,decreasing=TRUE),]
  outFile <- sprintf("%s/DMRs_cCRE_overlap_counts_%s.txt",outDir,dt)
  write.table(y,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  p <- ggplot(data = regions, aes(x=ENCODE.cCRE.Type)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("cis-Regulatory Element")
  pdfFile <- sprintf("%s/DMRs_cCRE_overlap_count_%s.pdf",
                     dmrDir, dt)
  ggsave(p,file=pdfFile)
  
  ## complete TE
  # TE class
  x <- table(TEregions$TE.repClass)
  y <- data.frame(
    TE.repClass=names(x),
    numDMRs_overlapping=as.numeric(x),
    pctOverlap=(as.numeric(x)/nrow(TEregions))
  )
  y <- y[order(y$numDMRs_overlapping,decreasing=TRUE),]
  outFile <- sprintf("%s/DMRs_TEclass_overlap_counts_%s.txt",outDir,dt)
  write.table(y,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  p <- ggplot(data = TEregions, aes(x=TE.repClass)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("cis-Regulatory Element")
  pdfFile <- sprintf("%s/DMRs_TEclass_overlap_count_%s.pdf",
                     dmrDir, dt)
  ggsave(p,file=pdfFile)
  
  # TE family
  x <- table(TEregions$TE.repFamily)
  y <- data.frame(
    TE.repFamily=names(x),
    numDMRs_overlapping=as.numeric(x),
    pctOverlap=(as.numeric(x)/nrow(TEregions))
  )
  y <- y[order(y$numDMRs_overlapping,decreasing=TRUE),]
  outFile <- sprintf("%s/DMRs_TEfamily_overlap_counts_%s.txt",outDir,dt)
  write.table(y,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  p <- ggplot(data = TEregions, aes(x=TE.repFamily)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("cis-Regulatory Element")
  pdfFile <- sprintf("%s/DMRs_TEfamily_overlap_count_%s.pdf",
                     dmrDir, dt)
  ggsave(p,file=pdfFile)
  
  ## simplified cCRE
  x <- table(simp_regions$ENCODE.cCRE.Type)
  y <- data.frame(
    ENCODE.cCRE.Category=names(x),
    numDMRs_overlapping=as.numeric(x),
    pctOverlap=(as.numeric(x)/nrow(simp_regions))
  )
  y <- y[order(y$numDMRs_overlapping,decreasing=TRUE),]
  outFile <- sprintf("%s/DMRs_cCRE_overlap_simplified_counts_%s.txt",outDir,dt)
  write.table(y,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  p <- ggplot(data = simp_regions, aes(x=ENCODE.cCRE.Type)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("cis-Regulatory Element")
  pdfFile <- sprintf("%s/DMRs_cCRE_overlap_simplified_count_%s.pdf",
                     dmrDir, dt)
  ggsave(p,file=pdfFile)
  
  ## simplified TE
  # TE class
  x <- table(simp_TEregions$TE.repClass)
  y <- data.frame(
    TE.repClass=names(x),
    numDMRs_overlapping=as.numeric(x),
    pctOverlap=(as.numeric(x)/nrow(simp_TEregions))
  )
  y <- y[order(y$numDMRs_overlapping,decreasing=TRUE),]
  outFile <- sprintf("%s/DMRs_cCRE_overlap_simplified_counts_%s.txt",outDir,dt)
  write.table(y,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  p <- ggplot(data = simp_TEregions, aes(x=TE.repClass)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("cis-Regulatory Element")
  pdfFile <- sprintf("%s/DMRs_TEclass_overlap_simplified_count_%s.pdf",
                     dmrDir, dt)
  ggsave(p,file=pdfFile)
  
  # TE family
  x <- table(simp_TEregions$TE.repFamily)
  y <- data.frame(
    TE.repFamily=names(x),
    numDMRs_overlapping=as.numeric(x),
    pctOverlap=(as.numeric(x)/nrow(simp_TEregions))
  )
  y <- y[order(y$numDMRs_overlapping,decreasing=TRUE),]
  outFile <- sprintf("%s/DMRs_TEfamily_overlap_simplified_counts_%s.txt",outDir,dt)
  write.table(y,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  keep <- unique(simp_TEregions$TE.repFamily[duplicated(simp_TEregions$TE.repFamily)]) # drop those only appeared once
  p <- ggplot(data = simp_TEregions[simp_TEregions$TE.repFamily %in% keep,], aes(x=TE.repFamily)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("cis-Regulatory Element (>1 overlap)") + 
    theme(axis.text.x = element_text(size = 5))
  pdfFile <- sprintf("%s/DMRs_TEfamily_overlap_simplified_count_%s.pdf",
                     dmrDir, dt)
  ggsave(p,file=pdfFile)
  
  
  ### GENCODE annotation ###
  message(sprintf("Reading in GENCODE annotation from: %s", geneFile))
  
  ## Load GENCODE db
  genes <- rtracklayer::readGFF(geneFile)
  
  ## GENCODE GRange Object - only protein coding
  genes_proteinCoding <- subset(genes, gene_type %in% "protein_coding" & type == "gene")
  genes_proteinCoding$TSS <- genes_proteinCoding$start
  genes_proteinCoding$TSS[which(genes_proteinCoding$strand=="-")] <- genes_proteinCoding$end[which(genes_proteinCoding$strand=="-")]
  
  geneGR_proteinCoding <- GRanges(genes_proteinCoding$seqid, 
                                  IRanges(genes_proteinCoding$TSS, genes_proteinCoding$TSS),
                                  name=genes_proteinCoding$gene_name
  ) 
  
  
  ## GENCODE GRange Object - all gene_type
  genes_all <- subset(genes, type == "gene")
  genes_all$TSS <- genes_all$start
  genes_all$TSS[which(genes_all$strand=="-")] <- genes_all$end[which(genes_all$strand=="-")]
  
  geneGR_all <- GRanges(genes_all$seqid, 
                        IRanges(genes_all$TSS, genes_all$TSS),
                        name=genes_all$gene_name
  ) 
  
  
  ## map PLS to nearest gene
  message("Map PLS to nearest genes...")
  pls <- grep("PLS",regions$ENCODE.cCRE.Type)
  pls <- regions[pls,]
  message(sprintf("\t%i PLS-overlapping DMRs",nrow(pls)))
  
  plsDMR <- GRanges(pls$ENCODE.cCRE.seqnames,
                    IRanges(pls$ENCODE.cCRE.start, pls$ENCODE.cCRE.end),
                    ID=pls$ENCODE.cCRE.ID)
  n <- nearest(plsDMR, geneGR_proteinCoding)
  pls$nearestTSS_proteinCoding <- geneGR_proteinCoding$name[n] 
  n <- nearest(plsDMR, geneGR_all)
  pls$nearestTSS_all <- geneGR_all$name[n] 
  
  outFile <- sprintf("%s/DMRs.ENCODE.PLS.nearestTSS_%s.txt",dmrDir, dt)
  write.table(pls,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  ## map ELS to nearest gene
  message("Map enhancers to nearest genes ...")
  els <- grep("E",regions$ENCODE.cCRE.Type)
  els <- regions[els,]
  message(sprintf("\t%i ELS-overlapping DMRs",nrow(els)))
  
  elsDMR <- GRanges(els$ENCODE.cCRE.seqnames,
                    IRanges(els$ENCODE.cCRE.start, els$ENCODE.cCRE.end),
                    ID=els$ENCODE.cCRE.ID)
  n <- nearest(elsDMR, geneGR_proteinCoding)
  els$nearestTSS_proteinCoding <- geneGR_proteinCoding$name[n] 
  n <- nearest(elsDMR, geneGR_all)
  els$nearestTSS_all <- geneGR_all$name[n] 
  
  outFile <- sprintf("%s/DMRs.ENCODE.ELS.nearestTSS_%s.txt",dmrDir, dt)
  write.table(els,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  
  ### Overlap DMR with Genes
  ol <- findOverlaps(
    ranges_dmrs, geneGR_all, minoverlap = 0
  )
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  
  message(sprintf("Found total %i OL; %i DMR X %i Genes", 
                  length(ol), uq_query, uq_sbj
  ))
  
  ## summarise overlap
  qry <- as.data.frame(ranges_dmrs[queryHits(ol)])
  colnames(qry) <- paste("DMR",colnames(qry),sep=".")
  qry_noHit <- as.data.frame(ranges_dmrs[-queryHits(ol)])
  
  sbj <- as.data.frame(geneGR_all[subjectHits(ol)])
  colnames(sbj) <- paste("GENE",colnames(sbj),sep=".")
  
  GENEregions <- cbind(qry,sbj)
  
  tmp <- cbind(qry_noHit, data.frame(matrix(nrow = nrow(qry_noHit), ncol = ncol(sbj))))
  colnames(tmp) <- colnames(GENEregions)
  tmp$GENE.name <- "noMatch"
  
  GENEregions <- rbind(GENEregions, tmp)
  
  simp_GENEregions <- GENEregions %>% 
    group_by(DMR.seqnames, DMR.start, DMR.end, DMR.width, DMR.strand, DMR.nCG, DMR.areaStat) %>% 
    summarise(
      GENE.name = toString(sort(unique(GENE.name))), 
      GENE.N_overlap = ifelse(GENE.name == "noMatch", 0, n()), 
    ) %>%
    arrange(.by_group = T)
  
  ## Save Gene overlap regions and simplified regions
  write.table(GENEregions,
              file=sprintf("%s/DMRs_Gene_overlap_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  
  write.table(simp_GENEregions,
              file=sprintf("%s/DMRs_Gene_overlap_simplified_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  rm(ol,qry,sbj)
  
  
  ### Overlap DMR with Aldinger CB CRE
  ## Load feCB_CRE db
  distal_enhancer_CB <- read.table(sprintf("%s/distalEnhancerLikeElements.bed", aldingerCB_CRE_dir), stringsAsFactors = F, header = F)
  proximal_enhancer_CB <- read.table(sprintf("%s/proximalEnhancerLikeElements.bed", aldingerCB_CRE_dir), stringsAsFactors = F, header = F)
  promoter_CB <- read.table(sprintf("%s/canonicalPromoterLikeElements.bed", aldingerCB_CRE_dir), stringsAsFactors = F, header = F)
  SE_CB <- read.table(sprintf("%s/superEnhancerLikeElements.bed", aldingerCB_CRE_dir), stringsAsFactors = F, header = F)
  
  distal_enhancer_CB$label <- "feCB_dELS"
  proximal_enhancer_CB$label <- "feCB_pELS"
  promoter_CB$label <- "feCB_PLS"
  SE_CB$label <- "feCB_SE"
  
  feCB_CRE <- rbind(distal_enhancer_CB, proximal_enhancer_CB, promoter_CB, SE_CB)
  print(table(feCB_CRE$label))
  
  feCB_CRE_Grange <- GRanges(feCB_CRE$V1,
                             IRanges((feCB_CRE$V2 + 1), feCB_CRE$V3) #convert bed 0-based to 1-based
  )
  mcols(feCB_CRE_Grange)$type <- feCB_CRE$label
  
  ## overlap with feCB_CRE
  ol <- findOverlaps(
    ranges_dmrs, feCB_CRE_Grange, minoverlap = 0
  )
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  
  message(sprintf("Found total %i OL; %i DMR X %i feCB_CREs", 
                  length(ol), uq_query, uq_sbj
  ))
  
  ## summarise
  qry <- as.data.frame(ranges_dmrs[queryHits(ol)])
  colnames(qry) <- paste("DMR",colnames(qry),sep=".")
  qry_noHit <- as.data.frame(ranges_dmrs[-queryHits(ol)])
  
  sbj <- as.data.frame(feCB_CRE_Grange[subjectHits(ol)])
  colnames(sbj) <- paste("feCB_CRE",colnames(sbj),sep=".")
  
  feCBregions <- cbind(qry,sbj)
  
  tmp <- cbind(qry_noHit, data.frame(matrix(nrow = nrow(qry_noHit), ncol = ncol(sbj))))
  colnames(tmp) <- colnames(feCBregions)
  tmp$feCB_CRE.type <- "noMatch"
  
  feCBregions <- rbind(feCBregions, tmp)
  
  simp_feCBregions <- feCBregions %>% 
    group_by(DMR.seqnames, DMR.start, DMR.end, DMR.width, DMR.strand, DMR.nCG, DMR.areaStat) %>% 
    summarise(
      feCB_CRE.type = toString(sort(unique(feCB_CRE.type))), 
      feCB_CRE.N_overlap = ifelse(feCB_CRE.type == "noMatch", 0, n())) %>%
    arrange(.by_group = T)
  
  ## Save feCB CRE overlap regions and simplified regions
  write.table(feCBregions,
              file=sprintf("%s/DMRs_feCB_overlap_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  
  write.table(simp_feCBregions,
              file=sprintf("%s/DMRs_feCB_overlap_simplified_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  rm(ol,qry,sbj)
  
  
  ### combine all simplified regions to one table
  message("Plotting Sankey plot to show overlap distribution")
  ## Merge and save
  simp_merged <- simp_regions %>% 
    merge(simp_feCBregions) %>%
    merge(simp_GENEregions) %>% 
    merge(simp_TEregions)
  
  write.table(simp_merged,
              file=sprintf("%s/DMRs_combined_overlap_simplified_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  
  ## Plot to show label connections
  # set thresholds to reduce figure complexity
  TEclass_keep <- as.data.frame(table(simp_TEregions$TE.repClass)) %>% 
    filter(Freq > 30) %>%
    mutate(Var1 = as.character(Var1)) %>%
    pull(Var1)
  ENCODE.cCRE.Type_keep <- as.data.frame(table(simp_regions$ENCODE.cCRE.Type)) %>% 
    filter(Freq > 50) %>%
    mutate(Var1 = as.character(Var1)) %>%
    pull(Var1)
  
  # Fortmat data
  tmp <- simp_merged %>% 
    mutate(id = paste(DMR.seqnames, DMR.start, DMR.end, sep = "_"), DNAm = ifelse(DMR.areaStat > 0, "RL-SVZ_hypoDNAm", "RL-VZ_hypoDNAm")) %>%
    select(c("id", "DNAm", "ENCODE.cCRE.Type", "feCB_CRE.type", "GENE.name", "TE.repClass")) %>%
    mutate(ENCODE.cCRE.Type = ifelse(ENCODE.cCRE.Type %in% ENCODE.cCRE.Type_keep, paste("ENCODE.cCRE.Type", ENCODE.cCRE.Type, sep = "_"), "ENCODE.cCRE.Type_Others")) %>%
    mutate(feCB_CRE.type = paste("feCB_CRE.type", feCB_CRE.type, sep = "_")) %>%
    mutate(GENE.name = ifelse(GENE.name != "noMatch", "GENE.name.Gene", paste("GENE.name", GENE.name, sep = "_"))) %>%
    mutate(TE.repClass = ifelse(TE.repClass %in% TEclass_keep, paste("TE.repClass", TE.repClass, sep = "_"), "TE.repClass_Others"))
  
  links_count <- NULL
  links_count <- rbind(links_count, as.data.frame(table(tmp$TE.repClass, tmp$ENCODE.cCRE.Type)))
  links_count <- rbind(links_count, as.data.frame(table(tmp$ENCODE.cCRE.Type, tmp$feCB_CRE.type)))
  links_count <- rbind(links_count, as.data.frame(table(tmp$feCB_CRE.type, tmp$DNAm)))
  links_count <- rbind(links_count, as.data.frame(table(tmp$DNAm, tmp$GENE.name)))
  colnames(links_count) <- c("source", "target", "value")
  
  nodes <- data.frame(
    name=c(as.character(links_count$source), 
           as.character(links_count$target)) %>% unique()
  )
  
  links_count$IDsource <- match(links_count$source, nodes$name)-1 
  links_count$IDtarget <- match(links_count$target, nodes$name)-1
  
  # plot sankey
  sankey <- sankeyNetwork(Links = links_count, Nodes = nodes, Source = "IDsource", Target = "IDtarget", Value = "value", NodeID = "name", sinksRight=FALSE, nodeWidth=50, fontSize=16, nodePadding=0)
  htmlwidgets::saveWidget(sankey, 
                          sprintf("%s/DMRs_combined_overlap_simplified_SankeyPlot_%s.html",dmrDir,dt),
                          selfcontained = TRUE)
  
}


### main ###
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

















