# characterize DMRs that overlap Group 3 and 4 MB CNAs.
rm(list=ls())

library(readxl)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38.masked) # needed for genNullSeqs

source("../../utils_PaiLab.R")
source("../overlapEnrichment/getGRanges_OLenrichment.R")

suppFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/ms_SuppTable/RhombicLipEpigenome_SuppTables_260331.xlsx"
epPairs <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Sarropoulos2026/260409_DMR_overlap_EPs/Sarropoulos_RL_DMR_overlap_inferredEP_upregGenes_FIMOhits_annotated.txt"

dmrFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"
geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

negDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRoverlap2/260401/negs"
negDirEP <- "/home/rstudio/isilon/private/projects/FetalHindbrain/CNA_DMRs_Overlap/260413/negs_EP_DMRs"

numPerm <- 1000L

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/CNA_DMRs_Overlap"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
 dir.create(outDir)
}

logFile <- sprintf("%s/%s.log", outDir, dt)
tryCatch({
 sink(logFile, split=TRUE)


cat("Collecting Northcott amp/dels\n")
n2012 <- suppressMessages(getNorthcott2012_AmpsDels())
n2017 <- suppressMessages(getNorthcott2017_AmpsDels())

hg38 <- BSgenome.Hsapiens.UCSC.hg38.masked

 cat("Starting DMR CNA overlap analysis...\n")
df <- as.data.frame(read_excel(suppFile, sheet="TableS9"))

n2012$amps[["GISTIC_Amps_MB"]] <- NULL
n2012$dels[["GISTIC_Dels_MB"]] <- NULL
n2012$amps[["GISTIC_Amps_SHH"]] <- NULL
n2012$dels[["GISTIC_Dels_SHH"]] <- NULL

n2017$amps[["SHH_GISTIC_AMP"]] <- NULL
n2017$dels[["SHH_GISTIC_DEL"]] <- NULL

#### compare the overlap between the 2012 and 2017 amp/del lists.
#### let's start with the dels
###dels_2012 <- n2012$dels[["GISTIC_Dels_Group3"]]
###dels_2017 <- n2017$dels[["GRP3_GISTIC_DEL"]]
###ol <- findOverlaps(dels_2012, dels_2017)
###cat(sprintf("Group 3 MB Dels: %i in 2012, %i in 2017, %i 2012 overlap %i 2017\n", 
###length(dels_2012), length(dels_2017), 
###    length(unique(queryHits(ol))), 
###    length(unique(subjectHits(ol)))))

cat("Combining amp/dels from 2012 and 2017 Northcott papers...\n")
CNAs <- list()
CNAs$amps <- list()
CNAs$amps$G3 <- suppressWarnings(reduce(c(n2012$amps[["GISTIC_Amps-Group3"]], n2017$amps[["GRP3_GISTIC_AMP"]])))
cat(sprintf("Group 3 MB Amps: %i in 2012, %i in 2017, %i total after merging\n", 
    length(n2012$amps[["GISTIC_Amps-Group3"]]), 
    length(n2017$amps[["GRP3_GISTIC_AMP"]]), 
    length(CNAs$amps$G3)))
CNAs$amps$G4 <- suppressWarnings(reduce(c(n2012$amps[["GISTIC_Amps_Group4"]], n2017$amps[["GRP4_GISTIC_AMP"]])))
cat(sprintf("Group 4 MB Amps: %i in 2012, %i in 2017, %i total after merging\n", 
    length(n2012$amps[["GISTIC_Amps_Group4"]]), 
    length(n2017$amps[["GRP4_GISTIC_AMP"]]), 
    length(CNAs$amps$G4)))

CNAs$dels$G3 <- suppressWarnings(reduce(c(n2012$dels[["GISTIC_Dels_Group3"]], n2017$dels[["GRP3_GISTIC_DEL"]])))
cat(sprintf("Group 3 MB Dels: %i in 2012, %i in 2017, %i total after merging\n", 
    length(n2012$dels[["GISTIC_Dels_Group3"]]), 
    length(n2017$dels[["GRP3_GISTIC_DEL"]]), 
    length(CNAs$dels$G3)))
CNAs$dels$G4 <- suppressWarnings(reduce(c(n2012$dels[["GISTIC_Dels_Group4"]], n2017$dels[["GRP4_GISTIC_DEL"]])))
cat(sprintf("Group 4 MB Dels: %i in 2012, %i in 2017, %i total after merging\n", 
    length(n2012$dels[["GISTIC_Dels_Group4"]]), 
    length(n2017$dels[["GRP4_GISTIC_DEL"]]), 
    length(CNAs$dels$G4)))      

dmrs <- read.delim(dmrFile, header=T, stringsAsFactors = FALSE)
dmrs$dmr_name <- sprintf("%s-%s-%s", dmrs$chr, dmrs$start, dmrs$end)
dmrGR <- GRanges(dmrs[,1], IRanges(start=dmrs[,2], end=dmrs[,3]), dmr_name = dmrs$dmr_name)

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

cat("\n")
cat("-------------\n")
 cat("Count E-P DMRs overlapping CNAs and affected genes\n")
 cat("-------------\n")
 ep <- read.delim(epPairs, header=T, sep="\t")
 ep$dmr_name <- sprintf("%s-%s-%s", ep$DMR_chr, ep$DMR_start, ep$DMR_end)
 ep <- ep[!duplicated(ep$dmr_name),]
 ep_dmrGR <- GRanges(seqnames=ep$DMR_chr, ranges=IRanges(start=ep$DMR_start, end=ep$DMR_end))
 ep_dmrGR$dmr_name <- ep$dmr_name
 ep_dmrGR$gene <- ep$gene

megaGeneList <- c()
overlappingCNAs <- list()
numCores <- 10L

cat("Examining overlap of eDMRs with Group 3 and 4 MB CNAs and affected genes...\n")
negTable_allDMRs <- list()
negTable_EPDMRs <- list()
for (curMega in c("amps","dels")) {
    cnas <- CNAs[[curMega]]
    overlappingCNAs[[curMega]] <- list()
    cat("-------------\n")
    cat(sprintf("Group 3 and 4 MB %s:\n", curMega))
    cat("-------------\n")
    for (cur in names(cnas)){

        cat("\n---\nTesting all DMRs\n----\n")
        ol_sv <- getGRanges_OLenrichment(
            pos=dmrGR,tgtGR=cnas[[cur]], numPerm=numPerm, negDir=negDir,
            rngSeed=12345,genome=hg38, outDir=outDir,
            tgtName=paste("AllDMRs_CNA",curMega,cur, collapse="_"), verbose=TRUE
        )

        tmp <- quantile(ol_sv$overlap_negs, probs=c(0.25, 0.50, 0.75))
        negTable_allDMRs[[paste(curMega,cur,sep="_")]] <- cbind(length(cnas[[cur]]), ol_sv$overlap_pos, 
            tmp[1], tmp[2], tmp[3], ol_sv$pval)


        cat("\n *** \nTesting E-P linked DMRs\n *** \n")
        ol_sv <- getGRanges_OLenrichment(
            pos=ep_dmrGR,tgtGR=cnas[[cur]], numPerm=100, negDir=negDirEP,
            rngSeed=12345,genome=hg38, outDir=outDir,
            tgtName=paste("CNA",curMega,cur, collapse="_"), verbose=TRUE
        )
        tmp <- quantile(ol_sv$overlap_negs, probs=c(0.25, 0.50, 0.75))
        negTable_EPDMRs[[paste(curMega,cur,sep="_")]] <- cbind(length(cnas[[cur]]), ol_sv$overlap_pos, 
            tmp[1], tmp[2], tmp[3], ol_sv$pval)

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
colnames(cnadf)[8:12] <- paste0("CNA_", colnames(cnadf)[8:12])
write.table(cnadf, file=sprintf("%s/eDMR_CNA_overlapping_genes.txt", outDir), sep="\t", quote=F, row.names=F)

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

negTable_allDMRs_df <- do.call(rbind, lapply(names(negTable_allDMRs), function(x) {
    data.frame(CNA=x, numGR=negTable_allDMRs[[x]][1], overlap_pos=negTable_allDMRs[[x]][2], 
    overlap_neg_Q25=negTable_allDMRs[[x]][3], overlap_neg_Q50=negTable_allDMRs[[x]][4], overlap_neg_Q75=negTable_allDMRs[[x]][5],
    pvalue=negTable_allDMRs[[x]][6],
        test="AllDMRs")
}))

negTable_EPDMRs_df <- do.call(rbind, lapply(names(negTable_EPDMRs), function(x) {
    data.frame(CNA=x, numGR=negTable_EPDMRs[[x]][1], overlap_pos=negTable_EPDMRs[[x]][2], 
        overlap_neg_Q25=negTable_EPDMRs[[x]][3], overlap_neg_Q50=negTable_EPDMRs[[x]][4], overlap_neg_Q75=negTable_EPDMRs[[x]][5],
        pvalue=negTable_EPDMRs[[x]][6],
        test="EP_DMRs")
}))

both <- rbind(negTable_allDMRs_df, negTable_EPDMRs_df)

# plot a two-panel barplot of the observed overlap, one for all DMRs and one for E-P linked DMRs, with error bars showing the median and interquartile range of the negative controls.
both$CNA <- factor(both$CNA, levels=unique(both$CNA))
# put three asterisk above the points with p <= 1e-3
p <- ggplot(both, aes(x=CNA, y=overlap_neg_Q50, fill=test)) +
 geom_bar(stat="identity", position=position_dodge()) +
 geom_errorbar(aes(ymin=overlap_neg_Q25, ymax=overlap_neg_Q75), position=position_dodge(0.9), width=0.25) +
 geom_point(aes(y=overlap_pos), position=position_dodge(0.9), color="red", size=3) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 labs(title="Overlap of DMRs with Group 3 and 4 MB CNAs", x="CNA", y="% DMRs overlapping CNA") +
 scale_fill_manual(values=c("AllDMRs"="grey50", "EP_DMRs"="grey80")) +
 theme(axis.text=element_text(size=18))
 p<- p + geom_text(data=both, aes(x=CNA, y=overlap_pos + 1, label=ifelse(pvalue <= 1e-3, "***", "")), position=position_dodge(0.9), size=5)
ggsave(p, file=sprintf("%s/DMR_CNA_overlap_barplot_withNegatives.pdf", outDir), width=10, height=6)
write.table(both, file=sprintf("%s/DMR_CNA_overlap_stats.txt", outDir), sep="\t", quote=F, row.names=F)


# genes overlapping DMR CNAs
cat(sprintf("Checking overlap of all DMRs with Group 3 and 4 MB CNAs and affected genes...\n"))
genes <- rtracklayer::import(geneFile)
genes <- subset(genes, gene_type == "protein_coding" & type == "gene")
geneTable <- list()
cnaDMR <- list()
dmrDirection <- list()
uq_genes <- NULL
for (curMega in names(CNAs)) {
    geneTable[[curMega]] <- list()
    cnaDMR[[curMega]] <- list()
    dmrDirection[[curMega]] <- list()
    for (cur in names(CNAs[[curMega]])) {
        dmrs_in_cna <- subsetByOverlaps(dmrGR, CNAs[[curMega]][[cur]])
        genes_in_cna <- findOverlaps(dmrs_in_cna, genes)
        geneTable[[curMega]][[cur]] <- cbind(
            as.data.frame(dmrs_in_cna[queryHits(genes_in_cna)]), 
            as.data.frame(genes[subjectHits(genes_in_cna)]))    

        geneTable[[curMega]][[cur]]$CNA <- paste(curMega, cur, sep="_")
        cat(sprintf("%s %s: %i gene overlaps\n", curMega, cur, nrow(geneTable[[curMega]][[cur]])))
        cat(sprintf("\t%i unique genes\n", length(unique(geneTable[[curMega]][[cur]]$gene_name))))
        uq_genes[[paste(curMega, cur, sep="_")]] <- length(unique(geneTable[[curMega]][[cur]]$gene_name))

        ol <- findOverlaps(dmrGR, CNAs[[curMega]][[cur]])
        cnaDMR[[curMega]][[cur]] <- as.data.frame(cbind(as.data.frame(dmrGR[queryHits(ol)]), 
            as.data.frame(CNAs[[curMega]][[cur]][subjectHits(ol)])))
        cnaDMR[[curMega]][[cur]]$CNA <- paste(curMega, cur, sep=" ")

        dmr_ol <- dmrs[which(dmrs$dmr_name %in% dmrGR$dmr_name[queryHits(ol)]),]
        dmr_ol <- dmr_ol[!duplicated(dmr_ol$dmr_name),]
        cat(sprintf("Num DMRs overlapping HARs: %i\n", length(unique(queryHits(ol)))))
        cat(sprintf("# hypoDMRs = %i, # hyperDMRs = %i\n", sum(dmr_ol$diff.Methy<0), sum(dmr_ol$diff.Methy>0)))
        dmrDirection[[curMega]][[cur]] <- cbind(sum(dmr_ol$diff.Methy<0), sum(dmr_ol$diff.Methy>0))
    }
}
uq_genes <- do.call(rbind, lapply(names(uq_genes), function(x) data.frame(CNA=x, gene=uq_genes[[x]])))
# make a barplot
uq_genes$CNA <- factor(uq_genes$CNA, levels=unique(uq_genes$CNA))
p <- ggplot(uq_genes, aes(x=CNA, y=gene, fill=CNA)) +
 geom_bar(stat="identity") +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 labs(title="Number of unique genes affected by DMRs overlapping Group 3 and 4 MB CNAs", x="CNA", y="Number of unique genes") +
 theme(axis.text=element_text(size=18)) +
 scale_fill_manual(values=c("amps_G3"="steelblue", "dels_G3"="steelblue4", "amps_G4"="coral", "dels_G4"="coral4"))
ggsave(p, file=sprintf("%s/DMR_CNA_overlap_genes_barplot.pdf", outDir), width=10, height=6)

browser()
cat("DMR direction for CNA overlap\n")
dmrDir <- do.call(rbind, lapply(names(dmrDirection), function(x) {
    do.call(rbind, lapply(names(dmrDirection[[x]]), function(y) {
        data.frame(CNA=paste(x,y,sep="_"), hypo=dmrDirection[[x]][[y]][1], hyper=dmrDirection[[x]][[y]][2])
    }))
}))
dmrDir$pctHyper <- (dmrDir$hyper / (dmrDir$hyper + dmrDir$hypo)) * 100
print(dmrDir)

# collapse into a table
geneTable_df <- do.call(rbind, lapply(geneTable, function(x) do.call(rbind, x)))
colnames(geneTable_df)[1:5] <- paste0("DMR_", colnames(geneTable_df)[1:5])
geneTable_df <- geneTable_df[,c(1:17, 33)]
geneTable_df <- geneTable_df[,-which(colnames(geneTable_df) %in% c("score", "phase", "gene_id", "gene_type"))]
colnames(geneTable_df)[6:ncol(geneTable_df)-2] <- paste0("gene_", colnames(geneTable_df)[6:ncol(geneTable_df)-2])
geneTable_df <- geneTable_df[!duplicated(geneTable_df),]
browser()
write.table(geneTable_df, file=sprintf("%s/allDMRs_overlapping_CNA_and_genes.txt", outDir), sep="\t", quote=F, row.names=F)

# collapse cnaDMR into a table
cnaDMR_df <- do.call(rbind, lapply(cnaDMR, function(x) do.call(rbind, x)))
colnames(cnaDMR_df)[1:5] <- paste0("DMR_", colnames(cnaDMR_df)[1:5])
colnames(cnaDMR_df)[6:10] <- paste0("CNA_", colnames(cnaDMR_df)[6:10])
write.table(cnaDMR_df, file=sprintf("%s/allDMRs_overlapping_CNA.txt", outDir), sep="\t", quote=F, row.names=F)

cat("Reading genes of interest\n")
neurodev_genes <- get_neurodevGenes()
g34_genes <- get_g34genes()  

 }, error=function(ex){
    print(ex)
}, finally = {
 sink()
})

