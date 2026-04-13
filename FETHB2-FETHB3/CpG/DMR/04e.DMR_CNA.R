# characterize DMRs that overlap Group 3 and 4 MB CNAs.
rm(list=ls())

library(readxl)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

suppFile <- "/Users/spai/Documents/FetalHindbrain_Epigenetics/RhombicLipEpigenome_SuppTables_260331.xlsx"
epPairs <- "/Users/spai/Documents/FetalHindbrain_Epigenetics/260409_Sarropoulos_DMR_overlap_EPs/Sarropoulos_RL_DMR_overlap_inferredEP_upregGenes_FIMOhits_annotated.txt"
###geneFile <- "/Users/spai/Documents/FetalHindbrain_Epigenetics/gencode.v44.basic.annotation.gtf.gz"

n2012_ampdel <- "/Users/spai/Documents/FetalHindbrain_Epigenetics/Northcott2012_AmpDels.Rdata"
n2017_ampdel <- "/Users/spai/Documents/FetalHindbrain_Epigenetics/Northcott2017_AmpDels.Rdata"

outDir <- "/Users/spai/Documents/FetalHindbrain_Epigenetics/CNA_DMRs_Overlap"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
 dir.create(outDir)
}

logFile <- sprintf("%s/%s.log", outDir, dt)
tryCatch({
 sink(logFile, split=TRUE)
 cat("Starting DMR CNA overlap analysis...\n")
df <- as.data.frame(read_excel(suppFile, sheet="TableS9"))
ep <- read.delim(epPairs, header=T, sep="\t",skip=1)

load(n2012_ampdel)
n2012_ampdel <- sv;

load(n2017_ampdel)
n2017_ampdel <- sv;
n2017_ampdel$amps$SHH_GISTIC_AMP <- NULL
n2017_ampdel$dels$SHH_GISTIC_DEL <- NULL


# count DMRs overlapping Group 3 and 4 MB CNAs.
# plot a bar graph of the number of DMRs overlapping each CNA, colored by whether they are in Group 3 or 4 MB.
cols <- colnames(df)[grepl("GISTIC", colnames(df))]
cnaOL <- list()
for (col in cols) {
#print(col)
 cnaOL[[col]] <- length(unique(df$name[df[,col]>0]))
}

cat(" Count DMRs overlapping CNAs:\n")
cnaOLdf <- data.frame(CNA=names(cnaOL), Count=unlist(cnaOL))
cnaOLdf$Group <- ifelse(grepl("G3", cnaOLdf$CNA), "Group 3", "Group 4")
p <- ggplot(cnaOLdf, aes(x=CNA, y=Count, fill=Group)) +
 geom_bar(stat="identity") +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 labs(title="Number of DMRs overlapping Group 3 and 4 MB CNAs", x="CNA", y="Number of DMRs")
 ggsave(p,file=sprintf("%s/DMR_CNA_overlap_barplot.png", outDir), width=8, height=6)

 cat("Count E-P DMRs overlapping CNAs and affected genes\n")
 ep$dmr_name <- sprintf("%s-%s-%s", ep$DMR_chr, ep$DMR_start, ep$DMR_end)
 ep <- ep[!duplicated(ep$dmr_name),]
 ep_dmrGR <- GRanges(seqnames=ep$DMR_chr, ranges=IRanges(start=ep$DMR_start, end=ep$DMR_end))
 ep_dmrGR$dmr_name <- ep$dmr_name
 ep_dmrGR$gene <- ep$gene

megaGeneList <- c()
overlappingCNAs <- list()
for (curMega in c("amps","dels")) {
    cnas <- n2017_ampdel[[curMega]]
    overlappingCNAs[[curMega]] <- list()
    cat("-------------\n")
    cat(sprintf("Group 3 and 4 MB %s:\n", curMega))
    cat("-------------\n")
    for (cur in names(cnas)){
        ol <- findOverlaps(ep_dmrGR, cnas[[cur]])
        cat(sprintf("%s: %i ranges, %i DMRs overlap\n", cur, length(cnas[[cur]]), length(unique(queryHits(ol)))))
        hit_dmrs <- ep_dmrGR$dmr_name[queryHits(ol)]

        ep_hits <- ep$gene[which(ep$dmr_name %in% hit_dmrs)]
        cat(sprintf("\t%i genes affected { %s }\n", length(unique(ep_hits)), paste(sort(unique(ep_hits)), collapse=", ")))
        megaGeneList <- c(megaGeneList, ep_hits)
        overlappingCNAs[[curMega]][[cur]] <- cbind(as.data.frame(ep_dmrGR[queryHits(ol)]), 
            as.data.frame(cnas[[cur]][subjectHits(ol)]))
        overlappingCNAs[[curMega]][[cur]]$gene <- ep_dmrGR$gene[queryHits(ol)]
    }
    cat("\n")
}
megaGeneList <- unique(megaGeneList)
cat(sprintf("Total unique genes affected by DMRs overlapping Group 3 and 4 MB CNAs: %i\n", length(megaGeneList)))
cat(sprintf("Genes: %s\n", paste(sort(megaGeneList), collapse=", ")))

cnadf <- do.call(rbind, lapply(overlappingCNAs, function(x) do.call(rbind, x)))
colnames(cnadf)[7]<- "predicted_target_gene"
colnames(cnadf)[8:13] <- paste0("CNA_", colnames(cnadf)[8:13])
write.table(cnadf, file=sprintf("%s/DMR_CNA_overlapping_genes.txt", outDir), sep="\t", quote=F, row.names=F)

cat("\n")
cat(sprintf("A total of %i DMRs overlap Group 3 and 4 MB CNAs, affecting %i unique genes.\n", 
    length(unique(cnadf$dmr_name)), 
    length(unique(cnadf$predicted_target_gene)))
)


#### checking gene coordinates
#### create a data.frame of TSS sites of genes
###gene.gr <- import(geneFile, format="gtf")
###gene.gr <- subset(gene.gr, gene_type == "protein_coding" & type == "transcript")
###tss.gr <- promoters(gene.gr, upstream=0, downstream=1)

browser()


 }, error=function(ex){
    print(ex)
}, finally = {
 sink()
})

