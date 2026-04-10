# diffEx RL VZ vs SVZ Aldinger RNA-seq
rm(list=ls())
require(edgeR)
require(ggplot2)
require(reshape2)

phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/input/Haldipur_RL_VZSVZ_counts/aldinger_rnaseq_0218_all.star_fc.metadata_rlvzsvz.txt"
datFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/input/Haldipur_RL_VZSVZ_counts/aldinger_rnaseq_0218_rl_0518.star_fc.counts_rlvzsvz.txt"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx"

# returns ggplot boxplot object
# gn = gene symbols
# 
plotGene <- function(gn, cpms, gp) {
    x <- melt(cpms[rownames(cpms) == gn,])
    x$group <- gp
    p <- ggplot(x, aes(y=value)) + geom_boxplot(aes(fill=gp))
    p <- p + ggtitle(gn) + ylim(0,14)
    p <- p + geom_hline(yintercept=0)
    p
}

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir)

pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
dat <- read.delim(datFile,sep="\t",h=T,as.is=T)
cat(sprintf("Counts, read %i genes, %i samples\n", 
    nrow(dat), ncol(dat)-1))
dat <- dat[!duplicated(dat$gene),]
cat(sprintf("After removing dups: %i genes\n", 
    nrow(dat)))
rownames(dat) <- dat$gene
dat <- dat[,-1]


# get sample name
x <- colnames(dat)
upos <- regexpr("_", x)
samp <- substr(x,2,upos[1]-1)
tis <- substr(x,upos[1]+1,nchar(x))
ID <- paste(samp,tis,sep="_")

message("matching pheno to data")
pheno$ID <- paste(pheno$donor, pheno$tissue,sep="_")
midx <- match(ID, pheno$ID)
if (all.equal(pheno$ID[midx], ID)!=TRUE){
    cat("mismatch")
    browser()
}
pheno <- pheno[midx,]

 idx <- which(pheno$age_pcw > 14)
    pheno <- pheno[idx, ]
    dat <- dat[, idx]
    cat(sprintf("After subsetting to age >14 pcw: %i samples\n", 
        ncol(dat)))

browser()

# edgeR
group <- factor(pheno$tissue)
y <- DGEList(counts=dat, group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
#y <- normLibSizes(y)
y <- calcNormFactors(y,method="TMM")
design <- model.matrix(~group)
message("estimating dispersion")
y <- estimateDisp(y,design)

message("performing gene-wise fit")
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
message("getting top tags")
tt <- lrt$table
tt$FDR <- p.adjust(tt$PValue,method="BH")
cat(sprintf("%i genes survive FDR correction\n", sum(tt$FDR < 0.05)))


write.table(tt,
	file=sprintf("%s/edgeR_RLVZvsSVZ_OlderThan14PCW_%s.txt",
	outDir,dt),sep="\t",col=T,row=T,,quote=F)

# plot BRCA1, WLS, OTX2, KAISO
pdfFile <- sprintf("%s/plotGenes.pdf", outDir)
pdf(pdfFile)
cpms <- cpm(y,log=TRUE)

x <-c("DCX","BRCA1","WLS","OTX2","ZBTB33","EOMES","CBFA2T2","CBFA2T3","SOX2","MKI67")
tryCatch({
for (gn in x){
    print(gn)
    p <- plotGene(gn, cpms, group)
    print(p)
}

library(ggrepel2)
# plot a volcano plot highlighting genes in x using ggrepel. Use ggplot2 for this. 
logFC <- tt$logFC
negLogPval <- -log10(tt$PValue)
df <- data.frame(logFC=logFC, negLogPval=negLogPval)
df$highlight <- ifelse(rownames(tt) %in% x, "yes",
"no")
p <- ggplot(df, aes(x=logFC, y=negLogPval))
p <- p + geom_point(aes(color=highlight))
p <- p + scale_color_manual(values=c("no"="grey","yes"="red"))
p <- p + geom_text_repel(data=subset(df, highlight=="yes"), 
    aes(label=rownames(subset(df, highlight=="yes"))), vjust=1.5)
p <- p + ggtitle("Volcano plot: RL vs VZ/SVZ diffEx")
ggsave(p,file=sprintf("%s/volcano_plot_RLVZvsSVZ_OlderThan14PCW_%s.png", outDir, dt), width=6, height=4, bg="white")



}, error=function(ex){
    print(ex)
}, finally={
    dev.off()
}
)