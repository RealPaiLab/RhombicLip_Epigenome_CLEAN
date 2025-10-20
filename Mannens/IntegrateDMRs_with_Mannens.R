# integrate DMRs with Mannens peak interaction data
rm(list=ls())
library(cicero)
library(GenomicRanges)

  dmrFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"

  mannensDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024_Cicero/250611"
  ciceroFile <- sprintf("%s/Mannens_cicero_conns_allchroms.qs",mannensDir)
  cdsFile <- sprintf("%s/Mannens_cicero_srat_cds_annotated_allchroms.qs",mannensDir)

  geneDef <- "/home/rstudio/isilon/src/gencode/GRCh38/gencode.v42.basic.annotation.gtf"

  outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/Mannens_Integration"
  dt <- format(Sys.Date(), "%y%m%d")
  outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/Mannens_Integration.log", outDir)
sink(logFile,split=TRUE)

tryCatch({

connsThresh <- 0.1

cat("* Reading DMRs from file\n")
cat("Input file: ", dmrFile, "\n")
dmrs <- read.delim(dmrFile, sep="\t",h=T)
cat(sprintf("Read %i DMRs from file\n", nrow(dmrs)))
dmr_GR <- makeGRangesFromDataFrame(dmrs, 
    seqnames.field = "chr", start.field = "start", end.field = "end",
    keep.extra.columns = TRUE)

cat("\n")
cat("* Reading Cicero connections from file\n")
cat("Cicero connections file: ", ciceroFile, "\n")
conns <- qs::qread(ciceroFile)
conns[,2] <- as.character(conns[,2])
cat(sprintf("Read %i Cicero connections from file\n", nrow(conns)))
cat(sprintf("Filtering Cicero connections with threshold %f\n", connsThresh))
conns <- subset(conns, coaccess >= connsThresh)
cat(sprintf("Found %i connections after filtering\n", nrow(conns)))

cat("\n")
cat("reading cell data set from file\n")
cds <- qs::qread(cdsFile)
  cat(sprintf("Read Cell Data Set with %i cells and %i genes\n", 
      ncol(cds), nrow(cds)))

###cat("* Reading annotated Cell Data Set object from file\n")
###cds <- qs::qread(cdsFile)
###cat(sprintf("Read Cell Data Set with %i cells and %i genes\n", 
###    ncol(cds), nrow(cds)))

cat("Finding DMR connections for each chromosome\n")
dmr_conns <- list()
t0 <- Sys.time()
for (chr in unique(seqnames(dmr_GR)) ){
  print(chr)
  cur <- subset(conns, 
      grepl(sprintf("%s-",chr), Peak1) | 
      grepl(sprintf("%s-",chr), Peak2)
  )
  cat(sprintf("\tFound %i connections for chromosome %s\n", 
      nrow(cur), chr))

  getOverlaps <- function(colName) {
    cur2 <- cur
    x <- strsplit(cur2[,colName], "-")
    cur2$chrom <- sapply(x, function(y) y[1])
    cur2$start <- as.numeric(sapply(x, function(y) y[2]))
    cur2$end <- as.numeric(sapply(x, function(y) y[3]))
  
    gr <- makeGRangesFromDataFrame(cur2,
    seqnames.field="chrom", start.field="start", end.field="end",
    keep.extra.columns = TRUE)

    overlaps <- findOverlaps(gr, dmr_GR)
    if (length(overlaps) == 0) {
      return(NULL)
    } else {
      dmr_overlaps <- dmr_GR[subjectHits(overlaps)]
      conns_overlaps <- gr[queryHits(overlaps)]
     
      combined <- data.frame(
        DMR_chr = seqnames(dmr_overlaps),
        DMR_start = start(dmr_overlaps),
        DMR_end = end(dmr_overlaps),
        Peak1 = conns_overlaps$Peak1,
        Peak2 = conns_overlaps$Peak2,
        coaccess = conns_overlaps$coaccess,
        stringsAsFactors = FALSE
      )
      return(combined)
    }
  }

  comb_pk1 <- getOverlaps("Peak1")
  comb_pk2 <- getOverlaps("Peak2")

  dmr_conns[[chr]] <- rbind(
    if (!is.null(comb_pk1)) comb_pk1 else data.frame(),
    if (!is.null(comb_pk2)) comb_pk2 else data.frame()
  )
}
print(Sys.time() - t0) # took 2.5 min on Pai Lab server

cat("Combining\n")
dmr_conns <- do.call(rbind, dmr_conns)
dmr_conns <- dmr_conns[!duplicated(dmr_conns),]
cat(sprintf("Found %i DMR connections in total\n", nrow(dmr_conns)))

cat("Writing DMR connections to file\n")
outFile <- sprintf("%s/Mannens_DMR_conns_%s.txt", outDir, dt)
write.table(dmr_conns, file=outFile, sep="\t", 
    quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("Finding unique DMRs from connections\n")
uq <- unique(dmr_conns[,c("DMR_chr", "DMR_start", "DMR_end")])
cat(sprintf("Found %i unique DMRs in total\n", nrow(uq)))
print(table(uq$DMR_chr))

cat("Finding overlaps of Cicero connections with TSS regions\n")
fd <- fData(cds)
fd$peak <- rownames(fd)

x <- merge(dmr_conns,fd,by.x=c("Peak1"), by.y=c("peak"), all.x=TRUE)
x2 <- merge(dmr_conns,fd,by.x=c("Peak2"), by.y=c("peak"), all.x=TRUE)
x <- rbind(x, x2)
x <- x[!is.na(x$gene),]
x <- x[!duplicated(x),]
cat(sprintf("Found %i DMR connections with TSS regions\n", nrow(x)))
y <- unique(x[,c("DMR_chr", "DMR_start", "DMR_end")])
cat(sprintf("Found %i unique DMRs with TSS regions\n", nrow(y)))
cat(sprintf("DMRs overlap with %i unique genes\n", 
    length(unique(x$gene))))

cat("Writing annotated DMR connections to file\n")
outFile <- sprintf("%s/Mannens_DMR_conns_annotated_%s.txt", outDir, dt)
write.table(x, file=outFile, sep="\t", 
    quote=FALSE, row.names=FALSE, col.names=TRUE)

},error=function(ex){
  print(ex)
},finally={
  cat("Closing log.\n")
  sessionInfo()
  sink(NULL)
})