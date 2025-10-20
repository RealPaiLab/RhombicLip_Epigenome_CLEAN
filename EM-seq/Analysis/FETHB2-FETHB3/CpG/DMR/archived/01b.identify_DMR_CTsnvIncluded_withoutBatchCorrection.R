# adapted from https://github.com/RealPaiLab/FetalHindbrain_Epigenetics/blob/master/FET_HB2/callDMR.R

rm(list=ls())

library(DSS)
library(ggplot2)

### dir config ###
rootDir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG"
cytoDir <- c("/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report", 
             "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/batch1_CpG_rerun")

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/DMRs/CTsnv_included/withoutBatchCorrection/%s",rootDir,dt)
logFile <- sprintf("%s/identify_DMR_CTsnvIncluded_withoutBatchCorrection_%s.log",outDir, dt)


main <- function() {
  message("DMR identification without batch correction")
  message(sprintf("Output will be written under: %s\n",outDir))
  
  ### analysis ###
  
  readBS <- function(inDir,grepPattern="",excludeSamp=c(), drop_hg38_alt = TRUE){
    fNames <- NULL
    inDirs <- NULL
    for (d in inDir) {
      message(sprintf("Searching files with %s pattern under %s ...", grepPattern, d))
      tmp <- dir(path=d,pattern=grepPattern)
      message(sprintf("%i files found", length(tmp)))
      fNames <- append(fNames, tmp)
      inDirs <- append(inDirs, rep(d, length(tmp)))
    }
    fNames <- fNames[grep(grepPattern,fNames)]
    sampNames <- fNames
    
    if (length(excludeSamp)>0){
      message(sprintf("Excluding {%s}",
                      paste(excludeSamp,collapse=",")))
      idx <- which(sampNames %in% excludeSamp)
      sampNames <- sampNames[-idx]
      fNames <- fNames[-idx]
      inDirs <- inDirs[-idx]
    }
    
    # test loci
    message(sprintf("about to read bismark from \n- %s", paste(paste(inDirs,fNames,sep="/"), collapse = "\n- ")))
    t0 <- Sys.time()
    bs <- bsseq::read.bismark(
      files=paste(inDirs,fNames,sep="/"),
      colData=DataFrame(inDir = inDirs, row.names = stringr::str_split(fNames, "\\.cytosine_report", simplify = T)[,1]),
      rmZeroCov=FALSE,
      strandCollapse=TRUE,
      verbose=TRUE 
    )
    t1 <- Sys.time()
    
    if (drop_hg38_alt) {
      bs <- bs[bs@rowRanges@seqnames %in% paste0("chr", c(as.character(1:22), "X", "Y"))] # remove spike-in and alt chrs
    } else {
      bs <- bs[bs@rowRanges@seqnames %in% c("pUC19", "J02459.1")] # only remove spike-in
    }
    
    return(list(files=fNames,bs=bs))
  }
  
  t0 <- Sys.time()
  message("Reading cytosine report files")
  bsObj <- readBS(inDir = cytoDir, 
                  grepPattern = "\\.txt\\.gz$")
  print(Sys.time()-t0)
  
  vzID <- c("FETHB2_0004_01_LB01-01_240206_A00469_0627_AHWVC5DSX7_1_TCTACGCA-GGCTATTG",
            "FETHB2_0005_01_LB01-01_240206_A00469_0627_AHWVC5DSX7_1_CTCAGAAG-AACTTGCC",
            "FETHB2_0002_01_LB01-01_230327_A00469_0456_BHHWCCDSX5_1_CATGAGCA-GCAAGATC",
            "FETHB2_0003_01_LB02-01_230327_A00469_0456_BHHWCCDSX5_1_AGGATAGC-AACCTTGG"
  )
  svzID <- c("FETHB2_0004_01_LB02-01_240206_A00469_0627_AHWVC5DSX7_1_GCAATTCC-TGTTCGAG",
             "FETHB2_0005_01_LB02-01_240206_A00469_0627_AHWVC5DSX7_1_GTCCTAAG-TGGTAGCT",
             "FETHB2_0001_01_LB01-01_230327_A00469_0456_BHHWCCDSX5_1_GAATCCGT-TGGAGTTG",
             "FETHB2_0002_01_LB02-01_230327_A00469_0456_BHHWCCDSX5_1_ATCTGACC-AAGTCGAG",
             "FETHB2_0003_01_LB01-01_230327_A00469_0456_BHHWCCDSX5_1_TCCTCATG-AGGTGTAC"
  )
  
  message("Performing DML tests")
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
                  minCG=4,
                  dis.merge=100,
                  pct.sig=0.5
  )
  
  dmrs_neg <- callDMR(dmlTest.sm,
                      delta=0,
                      p.threshold=1, 
                      minlen=50,
                      minCG=4,
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
  
  dmrs$status <- factor(ifelse(dmrs$diff.Methy > 0, 'Hypo', 'Hyper'), 
                        levels = c('Hyper', 'Hypo'))
  print(table(dmrs$status))
  
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
