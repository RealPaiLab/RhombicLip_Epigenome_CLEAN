# call DMRs

### NOTE: used oicr_hpc to generate the results. DSS_2.48.0 was installed with conda
### In HPC cmd:
### conda activate dss
### qsub -P pailab -V -cwd -b y -N dssdmr -M xsun@oicr.on.ca -m ea -l h_rt=1:0:0:0,h_vmem=10G -pe smp 32 Rscript --no-save ./01b.identify_DMR_CTsnvExcluded_withoutBatchCorrection.R

rm(list=ls())

library(DSS)
library(ggplot2)
require(reshape2)

source("../../utils_PaiLab.R")
source("getM_AllSamples_GRanges.R")

### dir config ###
rootDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG"
cytoDir <- get_configs("CPG_REPORT_DIR")
phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv"

dmrDir <- sprintf("%s/DMRs/CTsnv_excluded/withoutBatchCorrection", rootDir)

dmrFile <- sprintf("%s/250624/DMRs.csv", dmrDir)
#dmlFile <- sprintf("%s/250624/DMLs.Rdata", dmrDir)

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",dmrDir,dt)
logFile <- sprintf("%s/identify_DMR_CTsnvExcluded_withoutBatchCorrection_%s.log",
                   outDir, dt)

#dmrFile <- sprintf("%s/%s/DMRs.csv", outDir, dt)

  message("DMR identification without batch correction")
  message(sprintf("Output will be written under: %s\n",outDir))

  if (!dir.exists(outDir))  {
    dir.create(outDir, recursive = FALSE)
  }
  

sink(logFile,split=TRUE)
tryCatch({

  ### Analysis ###
  ## Load cytosine reports ##
  message("Reading cytosine report files")
  pheno <- read.delim(phenoFile, 
                        header = TRUE, 
                        stringsAsFactors = FALSE
  )
  dpos <- regexpr("-",pheno$library_ID)
  pheno$sample <- substr(pheno$library_ID, 1, dpos-1)

  if (file.exists(dmrFile)) {
    message(sprintf("DMR file already exists: %s", dmrFile))
    dmrs <- read.delim(dmrFile, sep="\t",header = TRUE, stringsAsFactors = FALSE)
    message(sprintf("Total DMRs: %d", nrow(dmrs)))
 
  } else {
    cat("Calling DMRs...\n")

    bsObj <- readBS(inDir = cytoDir, 
                  grepPattern = "\\.txt\\.gz$")
  
    browser()  
    ## Identify DMRs ##
    message("Performing DML tests")
    t0 <- Sys.time()

    dmlTest.sm <- DMLtest(bsObj$bs, 
                        group1=pheno$sample[which(pheno$ROI=="VZ")],
                        group2=pheno$sample[which(pheno$ROI=="SVZ")],
                        smoothing=TRUE, 
                        smoothing.span=500,
                        ncores=4
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
  
  pheno <- read.delim(phenoFile, 
                        header = TRUE, 
                        stringsAsFactors = FALSE
  )
  write.table(dmrs, file=sprintf("%s/DMRs.csv",outDir),
              sep="\t",col=T,row=F,quote=F)

  browser()
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
  ## Summary ##
  total_dmls <- nrow(dmlTest.sm)
  total_dmrs <- nrow(dmrs)

  message(sprintf("%d DMLs with p-value < 0.05 out of %d DMLs", 
                  sum(dmlTest.sm$pval < 0.05, na.rm = TRUE), 
                  total_dmls

                  )
          )
  message(sprintf("%d DMLs with q-value (FDR) < 0.05 out of %d DMLs", 
                  sum(dmlTest.sm$fdr < 0.05, na.rm = TRUE), 
                  total_dmls
                  )
          )
  message(sprintf("%d DMRs obtained out of %d DMLs", total_dmrs, total_dmls))
  
  ## Plot ##
  par(mfrow = c(1, 2))
  # p-value distribution of DMLs #
  p <- ggplot(dmlTest.sm, aes(x=pval)) +
    geom_bar(stat = 'bin', width = 0.1) +
    geom_vline(xintercept = 0.05, linetype = 'dashed') +
    xlab('Nominal p-value') +
    theme_minimal()
  outFile <- sprintf("%s/DML_pDistribution.png",outDir)
  ggsave(p, file=outFile)
  
  # DMR methylation status between tissues #
  dmrs$status <- factor(ifelse(dmrs$diff.Methy > 0, 'Hypo', 'Hyper'), 
                        levels = c('Hyper', 'Hypo'))
  print(table(dmrs$status))

cat("* Plotting diff.Methy histogram\n")
p <- ggplot(dmrs, aes(x = diff.Methy)) + # fill = status)) +
  geom_histogram(colour="white", fill="grey20",binwidth=0.05) +
  theme_classic(base_size = 14) +
  labs(x = 'DNAm Difference (%)',y='Number of DMRs')
ggsave(
  filename = sprintf("%s/DMR_histogram.pdf", outDir),
  plot = p,
  width = 8, height = 8
)
} # end else

  message("DMR identification completed")
  message(sprintf("Results saved to: %s", outDir))
  
cat("made it past the first one")
# sample 5 DMRs with high diff.Methy
gr <- makeGRangesFromDataFrame(dmrs, 
                            keep.extra.columns = TRUE, 
                            seqnames.field = "chr", 
                            start.field = "start", 
                            end.field = "end",
                            
)
names(gr) <- paste(dmrs$chr, 
                      dmrs$start, 
                      dmrs$end, 
                      sep = "-")

df <- as.data.frame(gr)
df$name <- rownames(df)
id <- df[,c("name","diff.Methy")]
id <- id[!duplicated(id),]
id$full <- sprintf("%s\n(diffMethy = %.2f)", id$name, id$diff.Methy)


cat("* Plotting sample DMRs with positive diff.Methy\n")
idx <- sample(which(gr$diff.Methy > 0.3), 5, replace = FALSE) 
a <- getM_AllSamples_GRanges(
  cytoDir, gr[idx])
pcta <- as.data.frame(a$pctM)
dpos <- regexpr("-",rownames(pcta))
rownames(pcta) <- substr(rownames(pcta), 1, dpos-1)
pcta$sample <- rownames(pcta)
pcta <- merge(pcta, pheno[,c("sample","ROI")], 
  by="sample", all.x=TRUE)

# plot boxplot of sample DMRs
x <- melt(pcta, id.vars = c("sample", "ROI"))
x <- na.omit(x)
x$ROI <- factor(x$ROI, levels=c("VZ","SVZ"))
x <- merge(x, id, by.x = "variable", by.y = "name", all.x = TRUE)
p <- ggplot(x, aes(x = ROI, y = value, fill = ROI)) +
  geom_boxplot(width=0.3, outliers=FALSE) + 
  geom_jitter(alpha = 0.5, width=0.3) +
  facet_grid(~ full, scales = "free_y") +
  labs(title = "Boxplot of DMRs across samples",
       x = "Sample",
       y = "Percent Methylation")  +
  scale_fill_manual(values= c("VZ" = "#af8dc3", "SVZ" = "#7fbf7b")) + 
  theme_minimal(base_size = 18)

ggsave(
  filename = sprintf("%s/DMR_posDiffMethy.pdf", outDir),
  plot = p,
  width = 16, height = 4
)  

cat("* Plotting sample DMRs with negative diff.Methy\n")
idx <- sample(which(gr$diff.Methy < -0.3), 10, replace = FALSE) 
a <- getM_AllSamples_GRanges(
  cytoDir, gr[idx])
pcta <- as.data.frame(a$pctM)
dpos <- regexpr("-",rownames(pcta))
rownames(pcta) <- substr(rownames(pcta), 1, dpos-1)
pcta$sample <- rownames(pcta)
pcta <- merge(pcta, pheno[,c("sample","ROI")], by="sample", all.x=TRUE)

# plot boxplot of sample DMRs
x <- melt(pcta, id.vars = c("sample", "ROI"))
x <- na.omit(x)

x$ROI <- factor(x$ROI, levels=c("VZ","SVZ"))
x <- merge(x, id, by.x = "variable", by.y = "name", all.x = TRUE)

for (cur in unique(x$variable)) {
  cur2 <- x[which(x$variable == cur), ]
  print(cur)

p <- ggplot(cur2, aes(x = ROI, y = value*100, fill = ROI)) +
  geom_boxplot(width=0.1) + 
  geom_point() +
  labs(title = cur,
       x = "Sample",
       y = "Percent Methylation")  +
  scale_fill_manual(values= c("VZ" = "#af8dc3", "SVZ" = "#7fbf7b")) + 
  theme_minimal(base_size = 12) 

  

  ggsave(p,file=sprintf("%s.pdf",cur), 
         width = 4, height = 4)

}

cat("finished plotting sample DMRs with negative diff.Methy\n")
browser()

# SP you identified a discrepancy in the diff.Methy values between the DMRs and the pctM values. You were going to see if the match improves if you first compute numerator and denominator separately, then compute the diff.Methy values.

# SP you were also going to load the DML result and look at base-resolution methylation computation by DSS. See how the %M is computed and how that matches
 # what you get from getM_AllSamples_GRanges()

ggsave(
  filename = sprintf("%s/DMR_negDiffMethy.pdf", outDir),
  plot = p,
  width = 16, height = 4
)  

# plot diff.Methy sorted by chromosome and start position, and change colours of diffmethy by chromosome
dmrs$chr <- factor(dmrs$chr, levels = paste0("chr", c(1:22, "X", "Y")))
cols <- rainbow(length(levels(dmrs$chr)))
names(cols) <- levels(dmrs$chr)
dmrs <- dmrs[order(dmrs$chr, dmrs$start), ]
dmrs$spos <- paste(dmrs$chr, dmrs$start, sep = ":")
p <- ggplot(dmrs, aes(x = spos, y = diff.Methy, colour = chr,alpha=0.5)) +
  geom_point() +
  scale_colour_manual(values=cols) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Differential Methylation by Chromosome and Start Position",
       x = "Chromosome:Start Position",
       y = "Differential Methylation") +

  theme_minimal(base_size = 18) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(
  filename = sprintf("%s/DMR_diffMethy_byChr.pdf", outDir),
  plot = p,
  width = 16, height = 10
)

# for each chromosome bin the diff.methy values into 10Mbp bins and take the average. Then plot sorted by chromosome and start position, color coded by chromosomes
dmrs$bin <- floor(dmrs$start / 1000000) * 1000000
dmrs$bin <- paste(dmrs$chr, dmrs$bin, sep = ":")
dmrs$bin <- factor(dmrs$bin, levels = unique(dmrs$bin))
# take the average diff.Methy for each bin
dmrs_avg <- aggregate(diff.Methy ~ bin, data = dmrs, FUN = mean)
# now plot the average diff.Methy by bin, colour coded by chromosome
dmrs_avg$chr <- sub("chr", "", sapply(strsplit(as.character(dmrs_avg$bin), ":"), `[`, 1))
dmrs_avg$chr <- factor(dmrs_avg$chr, levels = c(1:22, "X", "Y"))
cols <- rep(c("grey70","grey40"), length(levels(dmrs_avg$chr)))
names(cols) <- levels(dmrs_avg$chr)
dmrs_avg <- dmrs_avg[order(dmrs_avg$chr), ]
# plot them as rectangles
p <- ggplot(dmrs_avg, aes(x = bin, y = diff.Methy, colour = chr, alpha=0.5)) +
  geom_bar(stat="identity") +
  scale_colour_manual(values=cols) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Average Differential Methylation by Chromosome and Bin",
       x = "Chromosome:Bin",
       y = "Avg % methylation increase in RL-SVZ\n(binned 1MB)")
# draw vertical lines at the start of each chromosome
for (i in 1:length(levels(dmrs_avg$chr))) {
  chr <- levels(dmrs_avg$chr)[i]
  chridx <- which(dmrs_avg$chr == chr)
  fpos <- chridx[1]
  lpos <- chridx[length(chridx)]
  p <- p + geom_vline(xintercept = fpos, 
                      linetype = "dashed", color = "grey80")
  p <- p + annotate("text", 
                      x = (fpos+lpos)/2, 
                      y = max(dmrs_avg$diff.Methy, na.rm = TRUE) * 0.9, 
                      label = chr,                                     
                      size = 5)                     

  p <- p + theme_minimal(base_size = 18) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    legend.position = "none")
ggsave(
  filename = sprintf("%s/DMR_avgDiffMethy_byChr.pdf", outDir),
  plot = p,
  width = 14, height = 6
)
}
}, error = function(e) {message(e)},
         warning = function(w) {message(w)},
         finally = {        
           sink(NULL)                    
         }
) 
