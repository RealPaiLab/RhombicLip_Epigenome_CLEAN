#### THIS CODE IS UNFINISHED: LOAD CPH IS NOT FUNCTIONING PROPERLY ###

rm(list=ls())

library(DSS)
library(ggplot2)

### dir config ###
rootDir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpH"
cytoDir <- c("/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/CpH")

dt <- format(Sys.Date(),"%y%m%d")
#outDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection/%s",rootDir,dt)
#logFile <- sprintf("%s/identify_DMR_CTsnvExcluded_withoutBatchCorrection_%s.log",outDir, dt)

### analysis ###

## Load data
readBS <- function(inDir,grepPattern="",splitPattern="\\.cytosine_report",excludeSamp=c(), drop_hg38_alt = TRUE){
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

  bs <- bsseq::read.bismark(
    files=paste(inDirs,fNames,sep="/"),
    colData=DataFrame(inDir = inDirs, row.names = stringr::str_split(fNames, splitPattern, simplify = T)[,1]),
    rmZeroCov=FALSE,
    strandCollapse=TRUE,
    verbose=TRUE
  )
  
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
                grepPattern = "\\.txt\\.gz$", splitPattern = "\\.CHN")
print(Sys.time()-t0)


tmp <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/CpH/FETHB2_0001_01_LB01-01_230327_A00469_0456_BHHWCCDSX5_1_GAATCCGT-TGGAGTTG.CHN.cytosine_report.txt.gz"
tmp <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/alignment/methyldackel/report/FETHB2_0004_01_LB01-01_240206_A00469_0627_AHWVC5DSX7_1_TCTACGCA-GGCTATTG.cytosine_report.txt.gz"
bs <- bsseq::read.bismark(
  files=tmp,
  colData=DataFrame(inDir = "test", row.names = "FETHB2_0001_01_LB01-01_230327_A00469_0456_BHHWCCDSX5_1_GAATCCGT-TGGAGTTG"),
  rmZeroCov =FALSE,
  strandCollapse = FALSE,
  verbose=TRUE
)

