# adapted from https://github.com/RealPaiLab/FetalHindbrain_Epigenetics/blob/master/FET_HB2/dmrAnnotate.R

## ---- Identification of ccREs which overlap with DMRs

rm(list=ls())

library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)


### config ###
dt <- format(Sys.Date(),"%y%m%d")

rootDir <- "/data/xsun/output/EMseq_FETHB3"
outDir <- sprintf("%s/DMRs/CTsnv_included/%s",rootDir, "240330")

ccre_dir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/GRCh38-cCREs.bed"
dmrDir <- "/data/xsun/output/EMseq_FETHB3/DMRs/CTsnv_included/240330"

logFile <- sprintf("%s/annotateDMR_%s.log",outDir, dt)

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
  
  
  ### analysis ###
  df_dmrs <- dmrs
  df_ccres <- ccres
  
  ranges_ccres <- GRanges(ccres$V1, IRanges(ccres$V2, ccres$V3))
  ranges_ccres$ID <- ccres$V5
  ranges_ccres$Type <- ccres$V6
  
  ranges_dmrs <- GRanges(df_dmrs$chr, 
                         IRanges(df_dmrs$start, df_dmrs$end))
  ranges_dmrs$diff.Methy <- df_dmrs$diff.Methy
  
  
  ol <- findOverlaps(
    ranges_dmrs, ranges_ccres, minoverlap = 0
  )
  
  uq_query <- length(unique(queryHits(ol)))
  uq_sbj <- length(unique(subjectHits(ol)))
  
  message(sprintf("Found total %i OL; %i DMR X %i cCRE", 
                  length(ol), uq_query, uq_sbj
  ))
  
  qry <- as.data.frame(ranges_dmrs[queryHits(ol)])
  colnames(qry) <- paste("DMR",colnames(qry),sep=".")
  sbj <- as.data.frame(ranges_ccres[subjectHits(ol)])
  colnames(sbj) <- paste("ENCODE.cCRE",colnames(sbj),sep=".")
  
  regions <- cbind(qry,sbj)
  write.table(regions,
              file=sprintf("%s/DMRs_cCRE_overlap_%s.txt",dmrDir,dt),
              sep="\t",col=T,row=F,quote=F)
  rm(ol,qry,sbj)
  
  
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
  
  
  
  ### load GENCODE annotation ###
  geneFile <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/gencode.v42.basic.annotation.gtf"
  message("Reading in GENCODE annotation from: %s", geneFile)
  
  genes <- rtracklayer::readGFF(geneFile)
  genes <- subset(genes, gene_type %in% "protein_coding" & type == "gene")
  
  genes$TSS <- genes$start
  genes$TSS[which(genes$strand=="-")] <- genes$end[which(genes$strand=="-")]
  
  geneGR <- GRanges(genes$seqid, 
                    IRanges(genes$TSS, genes$TSS),
                    name=genes$gene_name
  ) 
  
  
  ### pLS to nearest gene ###
  message("Map pLS to nearest genes...")
  pls <- grep("PLS",regions$ENCODE.cCRE.Type)
  pls <- regions[pls,]
  message(sprintf("\t%i PLS-overlapping DMRs",nrow(pls)))
  
  plsDMR <- GRanges(pls$ENCODE.cCRE.seqnames,
                    IRanges(pls$ENCODE.cCRE.start, pls$ENCODE.cCRE.end),
                    ID=pls$ENCODE.cCRE.ID)
  n <- nearest(plsDMR, geneGR)
  pls$nearestTSS <- geneGR$name[n] 
  outFile <- sprintf("%s/DMRs.ENCODE.PLS.nearestTSS_%s.txt",dmrDir, dt)
  write.table(pls,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  ### ELS to nearest gene ###
  message("Map enhancers to nearest genes ...")
  els <- grep("E",regions$ENCODE.cCRE.Type)
  els <- regions[els,]
  message(sprintf("\t%i ELS-overlapping DMRs",nrow(els)))
  
  elsDMR <- GRanges(els$ENCODE.cCRE.seqnames,
                    IRanges(els$ENCODE.cCRE.start, els$ENCODE.cCRE.end),
                    ID=els$ENCODE.cCRE.ID)
  n <- nearest(elsDMR, geneGR)
  els$nearestTSS <- geneGR$name[n] 
  outFile <- sprintf("%s/DMRs.ENCODE.ELS.nearestTSS_%s.txt",dmrDir, dt)
  write.table(els,file=outFile,sep="\t",col=T,row=F,quote=F)
  
  
  ### ELS to ABC gene ###
  #nasserFile <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/Nasser2021_ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.neuronal.txt.gz"
  
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


















