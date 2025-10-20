# adapted from https://github.com/RealPaiLab/FetalHindbrain_Epigenetics/blob/master/FET_HB2/callDMR.R

rm(list=ls())

library(DSS)
library(ggplot2)

### dir config ###
rootDir <- "/data/xsun/output/EMseq_FETHB3"
cytoDir <- "/data/xsun/EMseq/FETHB3"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/DMRs/CTsnv_excluded/%s",rootDir,dt)
message(sprintf("Output will be written under: %s\n",outDir))
if (! dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

### CT SNPs ###
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)

### analysis ###
readBS <- function(inDir,grepPattern="",excludeSamp=c()){
  fNames <- dir(path=inDir,pattern=grepPattern)
  fNames <- fNames[grep(grepPattern,fNames)]
  sampNames <- fNames
  
  if (length(excludeSamp)>0){
    message(sprintf("Excluding {%s}",
                    paste(excludeSamp,collapse=",")))
    idx <- which(sampNames %in% excludeSamp)
    sampNames <- sampNames[-idx]
    fNames <- fNames[-idx]
  }
  
  # test loci
  message("about to read bismark")
  t0 <- Sys.time()
  bs <- bsseq::read.bismark(
    files=paste(inDir,fNames,sep="/"),
    colData=NULL,
    rmZeroCov=FALSE,
    strandCollapse=TRUE,
    verbose=TRUE
  )
  t1 <- Sys.time()
  
  return(list(files=fNames,bs=bs))
}

### preprocessing
t0 <- Sys.time()
message("Reading cytosine report files")
bsObj <- readBS(inDir = cytoDir, 
                grepPattern = ".txt$")
print(Sys.time()-t0)

### remove C to T SNPs
# method1 
#message("Filtering out C to T SNPs.")
#dbsnp <- SNPlocs.Hsapiens.dbSNP155.GRCh38
#chromosomes <- c(as.character(1:22), "X", "Y")
#test <- bsObj$bs
#message(sprintf("%s: %i methylation loci before filtering C to T SNPs.", Sys.time(),length(test)))
#for (chromosome in chromosomes) {
#  gc(verbose = F)
#  snps <- snpsBySeqname(dbsnp, seqnames = chromosome)
#  seqlevelsStyle(snps) <- "UCSC"
#  ref_alt <- inferRefAndAltAlleles(snps, BSgenome.Hsapiens.UCSC.hg38)
#  is_c2t <- (ref_alt$ref_allele == "C") & sapply(ref_alt$alt_alleles, function(x) {"T" %in% x})
#  c2t_snps <- snps[is_c2t]
#  test <- test[! test %over% c2t_snps]
#  
#  message(sprintf("%s: %i methylation loci after filtering Chr-%s C to T SNPs.", Sys.time(), length(test), chromosome))
#}
#message(sprintf("%s: %i methylation loci after filtering all C to T SNPs.", Sys.time(), length(test)))
# method2
bsObj$bs <- bsObj$bs[bsObj$bs@rowRanges@seqnames %in% paste0("chr", chromosomes <- c(as.character(1:22), "X", "Y"))] # remove spike-in and alt chrs
dbsnp <- SNPlocs.Hsapiens.dbSNP155.GRCh38
tmp <- bsObj$bs@rowRanges
seqlevelsStyle(tmp) <- "NCBI"
snpsByOverlaps(dbsnp, tmp[tmp@seqnames == "1"]) -> t
seqlevelsStyle(t) <- "UCSC"
ref_alt <- inferRefAndAltAlleles(t, BSgenome.Hsapiens.UCSC.hg38)
c2t_snps <- t[sapply(ref_alt$alt_alleles, function(x) {"T" %in% x})]
#!!! TODO: there are just too many extremely low MAF SNPs. maybe use genomad to filter instead? continue on AME regardless this for now.Not worth the time.



### dml
vzID <- c("/data/xsun/EMseq/FETHB3/FETHB2_0004_01_LB01-01_240206_A00469_0627_AHWVC5DSX7_1_TCTACGCA-GGCTATTG.cytosine_report.txt",
          "/data/xsun/EMseq/FETHB3/FETHB2_0005_01_LB01-01_240206_A00469_0627_AHWVC5DSX7_1_CTCAGAAG-AACTTGCC.cytosine_report.txt"
          )
svzID <- c("/data/xsun/EMseq/FETHB3/FETHB2_0004_01_LB02-01_240206_A00469_0627_AHWVC5DSX7_1_GCAATTCC-TGTTCGAG.cytosine_report.txt",
           "/data/xsun/EMseq/FETHB3/FETHB2_0005_01_LB02-01_240206_A00469_0627_AHWVC5DSX7_1_GTCCTAAG-TGGTAGCT.cytosine_report.txt"
           )

t0 <- Sys.time()
dmlTest.sm <- DMLtest(bsObj$bs, 
                      group1=vzID, 
                      group2=svzID,
                      smoothing=TRUE, 
                      smoothing.span=500
                      )
print(Sys.time()-t0)
message("saving DML data")
save(dmlTest.sm, file=sprintf("%s/DMLs.Rdata",outDir))

dmrs <- callDMR(dmlTest.sm,
                delta=0,
                p.threshold=1e-05, 
                minlen=50,
                minCG=3,
                dis.merge=100,
                pct.sig=0.5
)

dmrs_neg <- callDMR(dmlTest.sm,
                    delta=0,
                    p.threshold=1, 
                    minlen=50,
                    minCG=3,
                    dis.merge=100,
                    pct.sig=0.5
)

write.table(dmrs, file=sprintf("%s/DMRs.csv",outDir),
            sep="\t",col=T,row=F,quote=F)
write.table(dmrs_neg, file=sprintf("%s/DMRs_background.csv",outDir),
            sep="\t",col=T,row=F,quote=F)

## ---- Visualize a DMR, echo=TRUE, message=FALSE, fig.width=8, fig.height=10----

#showOneDMR(dmrs[1,], bsObj$bs)


total_dmls <- nrow(dmlTest.sm)
total_dmrs <- nrow(dmrs)

paste(sum(dmlTest.sm$pval < 0.05, na.rm = TRUE), "DMLs with p-value < 0.05 ","out of ", total_dmls)
paste(sum(dmlTest.sm$fdr < 0.05, na.rm = TRUE), "DMLs with q-value (FDR) < 0.05 ", "out of ", total_dmls)
paste(total_dmrs, "DMRs obtained out of ", total_dmls)

par(mfrow = c(1, 2))

p <- ggplot(dmlTest.sm, aes(x=pval)) +
  geom_bar(stat = 'bin', width = 0.1) +
  geom_vline(xintercept = 0.05, linetype = 'dashed') +
  xlab('Nominal p-value') +
  theme_minimal()
outFile <- sprintf("%s/DML_volcano.png",outDir)
ggsave(p, file=outFile)

dmrs$status <- factor(ifelse(dmrs$diff.Methy > 0, 'Hyper', 'Hypo'), 
                      levels = c('Hyper', 'Hypo'))

p1 <- ggplot(dmrs, aes(x = diff.Methy, fill = status)) +
  geom_histogram() +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Methylation Status')

p2 <- ggplot(dmrs, aes(x = diff.Methy, y = status, fill = status)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2) +
  theme_minimal() +
  scale_color_brewer(palette = 'Dark2') +
  labs( x = 'Change in Differential Methylation', 
        y = 'Methylation Status', 
        color = 'Methylation Status'
  )

outFile <- sprintf("%s/DMR_violins.pdf",outDir)
pdf(outFile)
print(p1)
print(p2)
dev.off()
