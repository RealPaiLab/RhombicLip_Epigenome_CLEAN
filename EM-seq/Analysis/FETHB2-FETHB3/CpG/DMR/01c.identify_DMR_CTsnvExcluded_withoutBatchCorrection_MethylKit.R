### On HPC:
### conda activate MethylKit
### qsub -P pailab -V -cwd -b y -N mkdmr -M xsun@oicr.on.ca -m ea -l h_rt=1:0:0:0,h_vmem=10G -pe smp 32 Rscript --no-save 01c.identify_DMR_CTsnvExcluded_withoutBatchCorrection_MethylKit.R

rm(list=ls())

library(methylKit)
library(glue)

#setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")

inputDir <- get_configs("CPG_REPORT_DIR")

rootDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG")
outDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection/%s/MethylKitDMRs",
                  rootDir, 
                  get_configs("CPG_DMR_DATE"))
dt <- format(Sys.Date(),"%y%m%d")

logFile <- sprintf(
  "%s/identify_DMR_CTsnvExcluded_withoutBatchCorrection_MethylKit_%s.log",
  outDir, 
  dt
  )
  
main <- function() {
  message("Using MethylKit to identify DMRs (without batch effect correction)")
  message(sprintf("Output will be written under: %s\n",outDir))
  
  ### load cytosine reports ###
  files <- list.files(
    inputDir,
    pattern = "*.txt.gz$")
  
  sample.id <- extract_sampleId(files)
  
  message("Loading cytosine reports")
  myobj <- methRead(as.list(paste(inputDir, files, sep = "/")),
                    sample.id = as.list(sample.id),
                    assembly = "hg38",
                    treatment = get_TissueType(sample.id),
                    context = "CpG", 
                    pipeline = "bismarkCytosineReport",
                    mincov = 1 # use mincov 1 as CpGs will be grouped by windows
  )
  
  ### Tiling genome CpGs ###
  win.size <- 50 # used 50 as DSS has min DMR 50 bp length
  step.size <- 50
  message(sprintf("Merging genomic regions by windows (window size: %d; step: %d)", 
                  win.size, 
                  step.size
                  )
          )
  tiles <- tileMethylCounts(myobj, 
                            win.size=win.size, 
                            step.size=step.size, 
                            cov.bases = 4, 
                            mc.cores = 28
                            )
  message("Saving tiles")
  saveRDS(tiles, sprintf("%s/tiles_%d_%d_%s.RDS", outDir, win.size, step.size, dt))
  
  message("Filtering for windows appearing in all samples")
  meth <- unite(tiles, destrand=T)
  message(sprintf("%s windows preserved", nrow(meth)))
  
  ### DMR identification ###
  message("Testing for differentially methylated regions")
  diff <- calculateDiffMeth(meth, 
                            covariates = NULL,
                            overdispersion="MN",
                            test="Chisq",
                            mc.cores = 28
                            )
  
  message("Saving diff object ...")
  saveRDS(diff, sprintf("%s/diff_%d_%d_%s.RDS", outDir, win.size, step.size, dt))
  
  message("Extracting DMRs ...")
  delta <- 10
  q_thresh <- 0.05
  dmr <- getMethylDiff(diff, difference=delta, qvalue=q_thresh)
  print(head(dmr))
  message(sprintf("# %d DMRs identified.", nrow(dmr)))
  
  message("Saving DMRs ...")
  write.table(dmr, 
              sprintf("%s/dmr_%d_%d_%d_%s.tsv", 
                      outDir, 
                      win.size, 
                      step.size, 
                      delta, 
                      dt
                      ),
              col.names = T, row.names = F, 
              sep = "\t"
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

