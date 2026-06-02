# compares CpG density in DMRs and controls.
rm(list=ls())
library(BSgenome.Hsapiens.UCSC.hg38.masked) # needed for genNullSeqs
library(Biostrings)
library(tidyr)
library(IRanges)
library(doParallel)
library(ggplot2)

source("../../utils_PaiLab.R")
source("utils.R")

dmrFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"
negDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRoverlap2/260601/negs"
SAMPLE_FOR_PLOT <- 10L

set.seed(42) # for reproducibility

dt <- format(Sys.Date(),"%y%m%d")

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_CpGcompare"
if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)

logFile <- sprintf("%s/compareCGdensity.log",outDir)
sink(logFile, split=TRUE)
tryCatch({

negs <- dir(path=negDir, pattern="neg")
negs <- negs[grep("(?<!unfiltered)\\.bed", negs, perl = TRUE)] # to use only filtered ones
negs <- sprintf("%s/%s",negDir,negs)
cat(sprintf("Found %d negative files.\n", length(negs)))

posFile <- sprintf("%s/pos1.fa",negDir)
pos <- readDNAStringSet(posFile)
cat("counting in positive seqs\n")
pos_counts <- dinucleotideFrequency(pos)[, "CG"] / width(pos)
pos_counts <- cbind(pos_counts, "pos")

cat("Distribution of CpG count in positive sequences:\n")
print(summary(as.numeric(pos_counts[,1])))

# use Biostrings library to count number of CpGs in each sequence.
registerDoParallel(cores=10)
cat("counting in negative seqs\n")
neg_counts <- foreach(i=1:length(negs), .packages = c("rtracklayer","Biostrings")) %dopar% {
  if (i %% 100 == 0) cat(sprintf("Processing negative file %d/%d\n", i, length(negs)))
  negFile <- negs[i]
  negGR <- import(negFile)
  negSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, negGR)
  x <- dinucleotideFrequency(negSeq)[, "CG"] / width(negGR)
  x <- cbind(x, sprintf("neg%d", i))
  x
}

cat("\nDistribution of CpG count in negative sequences:\n")
neg_counts_df <- do.call(rbind, neg_counts)
print(summary(as.numeric(neg_counts_df[,1])))

plot_idx <- sample(1:length(neg_counts), min(SAMPLE_FOR_PLOT, length(neg_counts)), FALSE)
cat(sprintf("Sampling for CpG plotting {%s}\n", paste(plot_idx, collapse=", ")))

  mega <- do.call(rbind, c(list(pos_counts), neg_counts[plot_idx]))
  mega <- as.data.frame(mega)
  mega[,1] <- as.numeric(mega[,1])
  colnames(mega) <- c("CpG_count", "Group")

  # plot a violin and skinny boxplot of counts in positive and negative ones. 
  # Show one plot per negative.
  p <- ggplot(mega, aes(x=Group, y=CpG_count, fill=Group)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1, fill="white") +
    theme_bw() +
    theme(legend.position="none", base_size=18) +
    labs(title="CpG count in positive and negative sequences", x="Group", y="CpG count")
  suppressWarnings(ggsave(sprintf("%s/CpG_count_violin.pdf", outDir), plot=p, width=8, height=6))
 ##$# make a close-up version of the above showing ylim c(0,100)
 ##$ p <- p + ylim(0, quantile(mega[,1], probs=0.95)) +
 ##$   labs(title="CpG count in positive and negative sequences (zoomed)", x="Group", y="CpG count")
 ##$ suppressWarnings(ggsave(sprintf("%s/CpG_count_violin_zoomed.pdf", outDir), plot=p, width=8, height=6))
 # plot the same thing as a density plot. Plot positive as red and negative as grey.
  p_density <- ggplot(mega, aes(x=CpG_count, colour=Group)) +
    geom_density(alpha=0.5) +
    scale_colour_manual(values=c("pos"="red", "neg1"="grey", "neg2"="grey", "neg3"="grey", "neg4"="grey", "neg5"="grey", "neg6"="grey", "neg7"="grey", "neg8"="grey", "neg9"="grey", "neg10"="grey")) +
    theme_bw() +
    labs(title="Density of CpG counts in positive and negative sequences", x="CpG count", y="Density")
  suppressWarnings(ggsave(sprintf("%s/CpG_count_density.pdf", outDir), plot=p_density, width=8, height=6))


# for each negative set, match sequences to the positive set by deciles.
pos_data <- mega[mega$Group == "pos", ]
quantiles <- quantile(pos_data[,1], probs=seq(0,1,0.25))

cat("Quantiles for positive CpG counts:\n")
print(quantiles)

# paralellize the loop below to speed up the matching process. 
# For each negative set, find the indices of sequences that fall into the same quantiles as the positive set. 
# Then save those indices to a file and create a new list of matched negatives.
newNegs <- foreach (i=1:length(neg_counts)) %dopar% {
  neg_name <- sprintf("neg%d", i)
  curidx <- list()
  cur <- neg_counts[[i]]


  for (d in 1:(length(quantiles)-1)) {
    curidx[[d]] <- sample(which(cur >= quantiles[d] & cur < quantiles[d+1]), length(cur)/length(quantiles), replace=TRUE)
    #cat(sprintf("\t%d: %i - %i ; %i seqs\n", d, quantiles[d], quantiles[d+1], length(curidx[[d]])))
  }
  match_idx <- unique(unlist(curidx))
  if (max(match_idx)> length(cur)) stop(sprintf("Warning: index exceeds bounds for file %d\n", i))
  if (i %% 100 == 0) cat(sprintf("%i:\t%i pos -> %i neg match\n", i, length(pos), length(match_idx)))
  #cat("\n")
  
  #write.table(match_idx, file=sprintf("%s/matched_indices_neg%d.txt", outDir, i), row.names=FALSE, col.names=FALSE)
  
  tmp <- data.frame(CpG_count = as.numeric(cur[match_idx]), Group = sprintf("neg%d", i))
} 

browser()


pos_counts <- as.data.frame(pos_counts)
pos_counts[,1] <- as.numeric(pos_counts[,1])
colnames(pos_counts) <- c("CpG_count", "Group")

# now plot the cpg count for positive and the newNegs.
mega <- do.call(rbind, c(list(pos_counts), newNegs[plot_idx]))
mega <- as.data.frame(mega)

p <- ggplot(mega, aes(x=Group, y=CpG_count, fill=Group)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  theme_bw() +
  theme(legend.position="none", base_size=18) +
  labs(title="CpG count in positive and matched negative sequences", x="Group", y="CpG count")
ggsave(sprintf("%s/CpG_count_violin_matched.pdf", outDir), plot=p, width=8, height=6)

# now subsample the positive set to match the deciles of the overall negative set.
negmega <- do.call(rbind, newNegs)
negmega <- as.data.frame(negmega)
negmega[,1] <- as.numeric(negmega[,1])
neg_quantiles <- quantile(negmega[,1], probs=seq(0,1,0.25))
cat("Quantiles for overall negative CpG counts:\n")
print(neg_quantiles)

browser()
matched_idx <- list()
new_pos <- foreach (d=1:(length(neg_quantiles)-1)) %do% {
  idx <- which(pos_data[,1] >= neg_quantiles[d] & pos_data[,1] < neg_quantiles[d+1])
  if (length(idx) == 0) {
    cat(sprintf("Warning: no positive sequences in quantile %d\n", d))
    return(NULL)
  }
  matched_idx[[d]] <- idx
  pos_data[idx, ]
}

matched_idx <- unique(unlist(matched_idx))
cat(sprintf("Total matched positive sequences: %i\n", length(matched_idx)))
write.table(matched_idx, file=sprintf("%s/matched_indices_pos.txt", outDir), row.names=FALSE, col.names=FALSE)

new_pos <- do.call(rbind, new_pos)
new_pos$Group <- "pos_matched"

new_negs <- lapply(plot_idx,  function(x) {
    tmp <- as.data.frame(newNegs[[x]])
    tmp
})
new_negs <- do.call(rbind, new_negs)
colnames(new_negs)[1] <- "CpG_count"

# now plot the new_pos against the newNegs.
mega <- rbind(new_pos, new_negs)

p <- ggplot(mega, aes(x=Group, y=CpG_count, fill=Group)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  theme_bw() +  
  theme(legend.position="none", base_size=18) +
  labs(title="CpG count in matched positive and negative sequences", x="Group", y="CpG count")
ggsave(sprintf("%s/CpG_count_violin_pos_matched.pdf", outDir), plot=p, width=8, height=6)

p <- p + ylim(0, 100) +
  labs(title="CpG count in matched positive and negative sequences (zoomed)", x="Group", y="CpG count")
ggsave(sprintf("%s/CpG_count_violin_pos_matched_zoomed.pdf", outDir), plot=p, width=8, height=6)

}, error=function(e) {
  cat("Error in processing negative files:", conditionMessage(e), "\n")
  return(NULL)
}, finally={
  cat("Finished processing negative files.\n")
  sink()
}
)