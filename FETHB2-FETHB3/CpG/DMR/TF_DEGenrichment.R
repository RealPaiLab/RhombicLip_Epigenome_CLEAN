# are TFs binding DMRs more enriched in SVG DEGs than other TFs?
rm(list=ls())

hypoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"

hyperFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/ame.tsv"

bgFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/MOTIFS_Hendrikse2022_RL_activeGenes.txt"

rnaDEGfile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/231123/edgeR_RLVZvsSVZ_231123.txt"

hypo <- read.delim(hypoFile, sep="\t", h=T, as.is=T)
na_idx <- which(is.na(hypo$adj_p.value))
if (length(na_idx)>0){
    hypo <- hypo[-na_idx,]
}
hyper <- read.delim(hyperFile, sep="\t", h=T, as.is=T)
na_idx <- which(is.na(hyper$adj_p.value))
if (length(na_idx)>0){
    hyper <- hyper[-na_idx,]
}

hypo <- hypo$motif_ID
dpos <- regexpr("\\.",hypo)
hypo_tf <- substr(hypo,1,dpos-1)
hypo_tf[which(hypo_tf=="NDF1")] <- "NEUROD1"
hypo_tf[which(hypo_tf=="NDF2")] <- "NEUROD2"

bg <- read.delim(bgFile, sep="\t", h=F, as.is=T)[,1]
bg <- sub("MOTIF ","", bg)
dpos <- regexpr("\\.",bg)
bg <- substr(bg,1,dpos-1)
bg[which(bg=="NDF1")] <- "NEUROD1"
bg[which(bg=="NDF2")] <- "NEUROD2"
cat(sprintf("Found %i TFs binding hypo DMRs\n", length(hypo_tf)))
cat(sprintf("Background TFs: %i\n", length(bg)))

deg <- read.delim(rnaDEGfile, sep="\t", h=T, as.is=T)
deg_genes <- rownames(deg[which(deg$PValue < 0.05),])
other_genes <- setdiff(rownames(deg), deg_genes)
cat(sprintf("Found %i DEGs (p<0.05) between RL-VZ and RL-SVZ\n", length(deg_genes)))

# % hypoTFs that are DEGs
num_hypo_tf_in_deg <- length(intersect(deg_genes, hypo_tf))
pct_hypo_tf_in_deg <- num_hypo_tf_in_deg / length(hypo_tf) * 100
cat(sprintf("Hypo DMR TFs in DEGs: %i / %i (%.2f%%)\n", 
    num_hypo_tf_in_deg, length(hypo_tf), pct_hypo_tf_in_deg))

# % other TFs that are DEGs
other_tf <- setdiff(bg, hypo_tf)
num_other_tf_in_deg <- length(intersect(deg_genes, other_tf))
pct_other_tf_in_deg <- num_other_tf_in_deg / length(other_tf) * 100
cat(sprintf("Other TFs in DEGs: %i / %i (%.2f
%%)\n", 
    num_other_tf_in_deg, length(other_tf), pct_other_tf_in_deg))

# binomial test
# null hypothesis: pct_other_tf_in_deg
# alternative hypothesis: pct_hypo_tf_in_deg > pct_other_tf_in_deg
p_value <- binom.test(num_hypo_tf_in_deg, length(hypo_tf), 
    p=pct_other_tf_in_deg/100, alternative="greater")$p.value
cat(sprintf("Binomial test p-value: %.4e\n", p_value))  

# Fisher's exact test
#            in_DEG not_in_DEG
# hypo_TF      a        b
# other_TF     c        d
a <- num_hypo_tf_in_deg
b <- length(hypo_tf) - num_hypo_tf_in_deg
c <- num_other_tf_in_deg
d <- length(other_tf) - num_other_tf_in_deg
fisher_matrix <- matrix(c(a,b,c,d), nrow=2,
    dimnames=list(c("hypo_TF","other_TF"), c("in_DEG","not_in_DEG")))
fisher_result <- fisher.test(fisher_matrix, alternative="greater")
cat("Fisher's exact test results:\n")
print(fisher_result)
