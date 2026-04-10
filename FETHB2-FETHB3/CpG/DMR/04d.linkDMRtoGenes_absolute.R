# link absolute methylation of DMRs to gene expression
rm(list=ls())

library(ggplot2)
library(reshape2)
library(edgeR)

#CPG_DMR_FILE <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"

dmr2genes <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_link2Genes_ABC/251007/DMR_AnnotatedAll_251007.tsv"

phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/input/Haldipur_RL_VZSVZ_counts/aldinger_rnaseq_0218_all.star_fc.metadata_rlvzsvz.txt"
rnaFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/input/Haldipur_RL_VZSVZ_counts/aldinger_rnaseq_0218_rl_0518.star_fc.counts_rlvzsvz.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMRtoGenes_absolute"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!file.exists(outDir)) dir.create(outDir,recursive=FALSE)

logFile <- sprintf("%s/linkDMRtoGenes_absolute_%s.log", outDir, dt)
sink(logFile, split=TRUE)

tryCatch({
    dat <- read.delim(dmr2genes,sep="\t",header=TRUE,stringsAsFactors=FALSE)

    pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
    pheno <- unique(pheno$donor[which(pheno$age_pcw > 14)])
    cat(sprintf("Found %i unique donors > 14 pcw\n", length(pheno)))

    cat("Read gene expression and normalize\n")
    xpr <- read.delim(rnaFile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    upos <- regexpr("_",colnames(xpr[,-1]))
    samp <- substr(colnames(xpr[,-1]),2,upos[1]-1)
    xpr <- xpr[,c(1, which(samp %in% pheno)+1)]
    cat(sprintf("After subsetting to %i samples with age >13 pcw\n", ncol(xpr)-1))

    dge <- suppressWarnings(DGEList(counts=xpr[,-1], genes=xpr[,1]))
    keep <- filterByExpr(dge)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge, method="TMM")
    cpms <- cpm(dge, log=TRUE, prior.count=1)
    xpr <- data.frame(gene=dge$genes$genes, cpms, stringsAsFactors=FALSE)

    # get average expression in RL-VZ and RL-SVZ
    xpr_rlvz <- xpr[,grep("_RL_vz",colnames(xpr))]
    xpr_rlsvz <- xpr[,grep("_RL_svz",colnames(xpr))]
    xpr$avg_rlvz <- rowMeans(xpr_rlvz)
    xpr$avg_rlsvz <- rowMeans(xpr_rlsvz)


    # print a violin plot of dat$diff.Methy
    p <- ggplot(dat, aes(x="", y=diff.Methy*100)) + geom_violin(trim=FALSE, fill="lightgreen") +
        geom_boxplot(width=0.1, outlier.shape=NA) +
        ylim(-100,100) +
        ggtitle("Distribution of absolute methylation difference of DMRs") +
        ylab("Absolute methylation difference (%)") +
        xlab("")
    outFile <- sprintf("%s/DMR_absolute_methylation_diff_violin.png", outDir)
    ggsave(outFile, p, width=5, height=5)


    # plot absolute methylation vs gene expression
    merged <- merge(dat, xpr, by.x="ABC_gene", by.y="gene", all.x=FALSE, all.y=FALSE)
    # plot avg_rlvz as a function of meanMethy2. First bin meanMethy2 into quartiles .
    merged$meanMethy2_bin <- cut(merged$meanMethy2, breaks=quantile(merged$meanMethy2, probs=seq(0,1,0.25)), include.lowest=TRUE)
    p <- ggplot(merged, aes(x=meanMethy2_bin, y=avg_rlsvz)) +
        geom_boxplot(aes(fill=meanMethy2_bin))
    p <- p + ggtitle("Gene expression in RL-SVZ vs binned absolute methylation") +
        xlab("Mean absolute methylation in DMR (binned into quartiles)") +
        ylab("Average gene expression in RL-SVZ (log CPM)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    outFile <- sprintf("%s/DMR_absolute_methylation_vs_RL_SVZ_expression.png", outDir)
    ggsave(outFile, p, width=7, height=5)
    




    browser()
    
}, error = function(e) {
  message("Error in linking DMRs to genes:")
  message(e)
},finally ={
    sink()
})