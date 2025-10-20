# adapted from https://github.com/RealPaiLab/FetalHindbrain_Epigenetics/blob/master/FET_HB2/callDMR.R

rm(list=ls())

library(DSS)
library(ggplot2)

### dir config ###
rootDir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG"
cytoDir <- c("/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB2/output/alignment/methyldackel", 
             "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report")

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/DMRs/CTsnv_included/dssBatchCorrected/%s",rootDir,dt)
logFile <- sprintf("%s/identify_DMR_CTsnvIncluded_dssBatchCorrection_%s.log",outDir, dt)


main <- function() {
  message("DMR identification with batch correction by DSS::DMLfit.multiFactor")
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
  region_dict <- c(rep("VZ", length(vzID)), rep("SVZ", length(svzID)))
  region_dict <- factor(region_dict, levels = c("VZ", "SVZ"))
  names(region_dict) <- c(vzID, svzID)
  
  # design matrix for linear effects
  message("Preparing for desgin matrix")
  design <- as.data.frame(bsObj$bs@colData)
  design$inDir <- NULL
  design$region <- region_dict[rownames(design)] # should consider about pair relation? i.e. vz/svz from the same donor
  design$batch <- ifelse(grepl("^FETHB2_000[1-3]_", rownames(design)), "FirstBatch", "SecondBatch")
  print(model.matrix(~region+batch, design))
  
  # Fit linear model
  t0 <- Sys.time()
  DMLfit = DMLfit.multiFactor(bsObj$bs, design=design, formula=~region+batch)
  message(Sys.time()-t0)
  message("Saving DMLfit data")
  save(DMLfit, file=sprintf("%s/DMLfit.Rdata",outDir))
  
  # DMLtest correct for batch effects
  t0 <- Sys.time()
  message("Perfroming DMLtest SVZ vs VZ (adjusted for batch)")
  DMLtest.region = DMLtest.multiFactor(DMLfit, coef="regionSVZ") # test region effect
  print(Sys.time()-t0)
  message("Saving DMLtest.region data")
  save(DMLtest.region, file=sprintf("%s/DMLtest.region.Rdata",outDir))
  
  # callDMR
  dmrs <- callDMR(DMLtest.region,
                  p.threshold=1e-02,
                  minlen=50,
                  minCG=4,
                  dis.merge=100,
                  pct.sig=0.5
  )
  
  dmrs_neg <- callDMR(DMLtest.region,
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
  
  
  total_dmls <- nrow(DMLtest.region)
  total_dmrs <- nrow(dmrs)
  
  paste(sum(DMLtest.region$pvals < 0.05, na.rm = TRUE), "DMLs with p-value < 0.05 ","out of ", total_dmls)
  paste(sum(DMLtest.region$fdrs < 0.05, na.rm = TRUE), "DMLs with q-value (FDR) < 0.05 ", "out of ", total_dmls)
  paste(total_dmrs, "DMRs obtained out of ", total_dmls)
  
  par(mfrow = c(1, 2))
  
  p <- ggplot(DMLtest.region, aes(x=pvals)) +
    geom_bar(stat = 'bin', width = 0.1) +
    geom_vline(xintercept = 0.05, linetype = 'dashed') +
    xlab('Nominal p-value') +
    theme_minimal()
  outFile <- sprintf("%s/DML_nomialP.png",outDir)
  ggsave(p, file=outFile)
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
