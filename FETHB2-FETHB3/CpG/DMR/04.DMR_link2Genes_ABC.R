# filters DMRs for those overlapping H3K27ac peaks and links to genes using neuro ABC

source("../../utils_PaiLab.R")
require(GenomicRanges)
library(UpSetR)
library(ggplot2)

CPG_DMR_FILE <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv"
abcFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/Nasser-Neuronal-ABC_creTarget_hg38.bed"

mannensGeneLink <- "/home/rstudio/isilon/src/neurodev-genomics/multiome/Mannens_2024/downloaded_from_authors/Gene_links.bedpe"

predictedTFBS1 <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hyper_AME_activeInRL_HOCOMOCOv12_240712/sequences.tsv"
predictedTFBS2 <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/SVZ_diff_hypo_AME_activeInRL_HOCOMOCOv12_240712/sequences.tsv"

harGR <- getHARs()
#### write HARs to a bed file
###write.table(
###    data.frame(
###        chr = seqnames(harGR),
###        start = start(harGR),
###        end = end(harGR)
###    ), 
###    file = "hars_hg38.bed")

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_link2Genes_ABC"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/DMR_link2Genes_ABC_%s.log", outDir, dt)
sink(logFile, split=TRUE)
tryCatch({
cat("reading dmrs\n")
dmrs <- read.delim(
    CPG_DMR_FILE, 
    header = TRUE, sep="\t",
    stringsAsFactors = FALSE
)

dmrs$name <- paste(dmrs$chr, 
                      dmrs$start, 
                      dmrs$end, 
                      sep = "-"
)
dmrGR <- makeGRangesFromDataFrame(
    dmrs, 
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end"
)

abc <- read.delim(abcFile, header=FALSE, stringsAsFactors = FALSE)
abc <-  subset(abc, V1 %in% paste("chr", c(1:22, "X", "Y"),sep=""))
abcGR <- makeGRangesFromDataFrame(
    abc, 
    keep.extra.columns = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3"
)
# create a combined table of DMRs overlapping ABC
combined <- GenomicRanges::findOverlaps(dmrGR, abcGR)
x1 <- dmrs[queryHits(combined),]
x2 <- abc[subjectHits(combined),]
both <- cbind(x1, x2)

# print unique DMRs in overlap, and how many unique genes are linked
cat("\nDMRs overlapping ABC:\n")
cat(sprintf("Number of unique DMRs overlapping ABC: %d\n", length(unique(both$name))))
cat(sprintf("Number of unique genes linked to DMRs: %d\n", length(unique(both$V4))))

colnames(both)[which(colnames(both) == "V4")] <- "ABC_gene"

uqGenes <- unique(both$ABC_gene)
neurodev_genes <- get_neurodevGenes()
  both$IsNeurodev <- both$ABC_gene %in% neurodev_genes
  
cat("* Getting Group 3 & 4 MB associated genes\n")
g34_genes <- get_g34genes()  
both$IsG34gene <- both$ABC_gene %in% g34_genes  
  
both <- both[, -which(colnames(both) %in% c("V1","V2","V3"))]

# Now read in Mannens gene links
mannens <- read.delim(
    mannensGeneLink, 
    header = FALSE,
    stringsAsFactors = FALSE
)
mannensGR <- makeGRangesFromDataFrame(
    mannens, 
    keep.extra.columns = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3"
)
mannensGR <- resize(mannensGR, fix = "center", width = 1000)
# create a combined table of DMRs overlapping Mannens gene links
combined <- GenomicRanges::findOverlaps(dmrGR, mannensGR)
x1 <- dmrs[queryHits(combined),]
x2 <- mannens[subjectHits(combined),]
mannensBoth <- cbind(x1, x2)  
colnames(mannensBoth)[which(colnames(mannensBoth) == "V7")] <- "Mannens_gene"  
mannensBoth <- mannensBoth[!duplicated(mannensBoth), ]

cols <- c("chr", "start", "end", "length",
    "nCG","meanMethy1","meanMethy2","diff.Methy","areaStat")
x <- merge(both, mannensBoth, by = "name",all.x=TRUE, all.y=TRUE)
for (i in cols){
    na_idx <- which(is.na(x[,sprintf("%s.x", i)]))
    if (length(na_idx) == 0) next
    x[na_idx, sprintf("%s.x", i)] <- x[na_idx, sprintf("%s.y", i)]
 #   x <- x[, -which(colnames(x) %in% c(sprintf("%s.y", i)))]
}
colnames(x)<- sub("\\.x$", "", colnames(x))
idx <- grep("\\.y$", colnames(x))
x <- x[, -idx]
x <- x[,-which(colnames(x) %in% paste("V",c(1:6,8),sep=""))]
write.table(
    x, 
    file = sprintf("%s/Mannens_gene_links_%s.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)
cat("get fetal CB enhancers\n")
enh <- get_fetalCB_enh()
cat("filtering DMRs for enhancers\n")
# create a combined table of DMRs overlapping fetal CB enhancers
combined <- GenomicRanges::findOverlaps(dmrGR, enh)
olCB <- rep(FALSE, length(dmrGR))
olCB[queryHits(combined)] <- TRUE
dmrs$Overlaps_FetalCB_Enhancer <- olCB

both3 <- merge(x, dmrs, by = "name", all.x = TRUE, all.y = TRUE)
for (i in cols){
    print(i)
    na_idx <- which(is.na(both3[,sprintf("%s.x", i)]))
    if (length(na_idx) == 0) next
    both3[na_idx, sprintf("%s.x", i)] <- both3[na_idx, sprintf("%s.y", i)]
}
colnames(both3)<- sub("\\.x$", "", colnames(both3))
idx <- grep("\\.y$", colnames(both3))
both3 <- both3[, -idx]
colnames(both3)[which(colnames(both3) == "IsNeurodev")] <- "ABC_IsNeurodev"
colnames(both3)[which(colnames(both3) == "IsG34gene")] <- "ABC_IsG34gene"
both3 <- both3[!duplicated(both3),]

write.table(
    both3, 
    file = sprintf("%s/DMR_AnnotatedAll_%s.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)
 
sv <- getNorthcott2012_AmpsDels()
# make SIF file for Cytoscape
olList <- list()
for (x in c("amps", "dels")) {
    for (nm in names(sv[[x]])) {
        if (any(grep("SHH", nm))) next;
        olList[[nm]] <- get_ol_targets(dmrGR,sv[[x]][[nm]],
        nm)
        both3[, sprintf("Overlaps_N2012_%s", nm)] <- FALSE
        if (length(olList[[nm]]) == 0) next
        idx <- which(both3$name %in% olList[[nm]]$name)
        both3[idx, sprintf("Overlaps_N2012_%s", nm)] <- TRUE
    }
}

sv <- getNorthcott2017_AmpsDels()
# make SIF file for Cytoscape
olList <- list()
for (x in c("amps", "dels")) {
    for (nm in names(sv[[x]])) {
        if (any(grep("SHH", nm))) next;
        olList[[nm]] <- get_ol_targets(dmrGR,sv[[x]][[nm]],
        nm)
        both3[, sprintf("Overlaps_N2017_%s", nm)] <- FALSE
        if (length(olList[[nm]]) == 0) next
        
        idx <- which(both3$name %in% olList[[nm]]$name)
        both3[idx, sprintf("Overlaps_N2017_%s", nm)] <- TRUE
    }
}

write.table(
    both3, 
    file = sprintf("%s/DMR_AnnotatedAll_WithCNAs_%s.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)

olSV <- rowSums(both3[, grep("Overlaps_N2012|Overlaps_N2017", colnames(both3))]
, na.rm = TRUE) > 0
nm <- unique(both3$name[olSV])
has_SV <- unique(both3$name[which(olSV)])


cat("Final set of statistics\n")
cat(sprintf("Total unique DMRs: %d\n", 
length(unique(both3$name))))
idx <- which(!is.na(both3$ABC_gene))
cat(sprintf("Total unique DMRs overlapping ABC: %d\n", 
    length(unique(both3$name[idx]))))
idx <- which(!is.na(both3$Mannens_gene))
cat(sprintf("Total unique DMRs overlapping Mannens gene links: %d\n", 
    length(unique(both3$name[idx]))))
idx <- which(both3$Overlaps_FetalCB_Enhancer)
cat(sprintf("Total unique DMRs overlapping fetal CB enhancers: %d\n", 
    length(unique(both3$name[idx]))))
cat(sprintf("Total unique neurodev genes: %i\n", 
    length(unique(both3$ABC_gene[both3$ABC_IsNeurodev]))))
cat(sprintf("Total unique G34 MB genes: %i\n", 
    length(unique(both3$ABC_gene[both3$ABC_IsG34gene]))))
cat(sprintf("Total unique DMRs overlapping Northcott 2012 or 2017 SV: %d\n", 
    length(has_SV)))

write.table(
    both3, 
    file = sprintf("%s/DMR_AnnotatedAll_WithCNAs_%s.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)

# write SIF file.
x <- both3[,-c(2:10)]


nm <- x$name[which(x$ABC_IsNeurodev | x$ABC_IsG34gene)]
x <- x[which(x$name %in% nm),]
x <- x[!duplicated(x),]
tmp <- rep(FALSE, nrow(x))
tmp[which(x$name %in% has_SV)] <- TRUE
x$olSV <- tmp

cat("writing gene attributes\n")
geneAttrs <- x
geneAttrs <- geneAttrs[,c("ABC_gene","ABC_IsNeurodev",
    "ABC_IsG34gene","Mannens_gene")]
geneAttrs <- geneAttrs[!duplicated(geneAttrs),]
geneType <- rep(NA, nrow(geneAttrs))
geneType[which(geneAttrs$ABC_IsNeurodev)] <- "neurodev"
geneType[which(geneAttrs$ABC_IsG34gene)] <- "G34"
geneAttrs$geneType <- geneType
write.table(
    geneAttrs, 
    file = sprintf("%s/DMR_link2Genes_ABC_%s_geneAttrs.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)   

cat("writing DMR attributes\n")
dmrAttrs <- x
dmrAttrs <- dmrAttrs[,c("name","Overlaps_FetalCB_Enhancer","olSV")]
dmrAttrs <- dmrAttrs[!duplicated(dmrAttrs),]

write.table(
    dmrAttrs, 
    file = sprintf("%s/DMR_link2Genes_ABC_%s_dmrAttrs.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)

cat("writing SIF file\n")
sifOut <- x[,1:2]; 
sifOut <- sifOut[!duplicated(sifOut),]
sifOutFile <- sprintf("%s/DMR_link2Genes_ABC_%s.sif", outDir, dt)
if (file.exists(sifOutFile)) {
    file.remove(sifOutFile)
}
for (i in unique(sifOut[,2])){
    idx <- which(sifOut[,2] == i)
    if (length(idx) == 0) next
    cat(sprintf("%s\tregulates\t%s\n", i, 
        paste(sifOut[idx,1], collapse = "\t")),
        file = sifOutFile, append = TRUE)
}

cat("Now create links from genes to enhancers\n")
tfbs1 <- read.delim(predictedTFBS1, header=TRUE, stringsAsFactors = FALSE)
tfbs2 <- read.delim(predictedTFBS2, header=TRUE, stringsAsFactors = FALSE)
tfbs <- rbind(tfbs1, tfbs2)
tfbs <- subset(tfbs, class %in% "tp")[,c(2,4)]
cat(sprintf("Read %i predicted TFBS\n", nrow(tfbs)))
tfbs <- tfbs[!duplicated(tfbs),]
cat(sprintf("%i unique TF motifs\n", length(unique(tfbs$motif_ID))))
cat(sprintf("After removing duplicates, %i predicted TFBS\n", nrow(tfbs)))
tf <- unique(tfbs$motif_ID)
tfbs[,1] <- substr(tfbs[,1], 1, regexpr("\\.", tfbs[,1])-1)
cat(sprintf("After removing sub-motif info, %i unique TF motifs\n", length(unique(tfbs$motif_ID))))
tfbs <- tfbs[!duplicated(tfbs),]
cat(sprintf("After removing duplicate motif IDs, %i predicted TFBS\n", nrow(tfbs)))
tfbs[,2] <- sub(":", "-", tfbs[,2])
#split column 2 by "-"
tfbs_split <- do.call(rbind, strsplit(tfbs[,2], "-"))
tfbs_gr <- makeGRangesFromDataFrame(
    data.frame(
        chr = tfbs_split[,1],  
        start = as.integer(tfbs_split[,2]),
        end = as.integer(tfbs_split[,3]),
        motif_ID = tfbs[,1]
    ), keep.extra.columns = TRUE)

# now split sifOut[,2] by "-", make GRanges and find overlaps
sif_split <- do.call(rbind, strsplit(sifOut[,1], "-"))
sif_gr <- makeGRangesFromDataFrame(
    data.frame(
        chr = sif_split[,1],
        start = as.integer(sif_split[,2]),
        end = as.integer(sif_split[,3]),
        name = sifOut[,2]
    ), keep.extra.columns = TRUE)

combined <- GenomicRanges::findOverlaps(sif_gr, tfbs_gr)
x1 <- sifOut[queryHits(combined),]
x2 <- tfbs_gr[subjectHits(combined),]
tfbs_sif <- cbind(x1, as.data.frame(x2))
colnames(tfbs_sif)[1] <- c("enhancer")
tfbs_sif <- tfbs_sif[,c("ABC_gene","motif_ID")]
tfbs_sif <- tfbs_sif[!duplicated(tfbs_sif),]

TF2gene <- tfbs_sif[,c("motif_ID","ABC_gene")]

TF2gene <- subset(TF2gene,
    motif_ID %in% c(neurodev_genes, g34_genes))
cat(sprintf("After filtering for neurodev and G34 genes, %i links remain\n", nrow(TF2gene)))
write.table(
    TF2gene, 
    file = sprintf("%s/TFmotif_2_ABCgene_MotifIDIsNeurovOrG34_%s.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)

TF2gene <- subset(TF2gene,
    ABC_gene %in% c(neurodev_genes, g34_genes))
cat(sprintf("After filtering for ABCgene being neurodev or G34 genes, %i links remain\n", nrow(TF2gene)))
write.table(
    TF2gene, 
    file = sprintf("%s/TFmotif_2_ABCgene_srcTgt_IsNeurovOrG34_%s.tsv", outDir, dt), 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE
)

cat("create clean table for upset plot\n")
forup <- data.frame(
    name = dmrs$name,
    Overlaps_FetalCB_Enhancer = FALSE,
    ABC_gene = FALSE,
    ABC_IsNeurodev = FALSE,
    ABC_IsG34gene = FALSE,
    olSV = FALSE
)
forup[which(forup$name %in% both3$name[which(both3$Overlaps_FetalCB_Enhancer)]), 
    "Overlaps_FetalCB_Enhancer"] <- TRUE
forup[which(forup$name %in% both3$name[which(!is.na(both3$ABC_gene))]), 
    "ABC_gene"] <- TRUE
forup[which(forup$name %in% both3$name[which(both3$ABC_IsNeurodev)]), 
    "ABC_IsNeurodev"] <- TRUE
forup[which(forup$name %in% both3$name[which(both3$ABC_IsG34gene)]), 
    "ABC_IsG34gene"] <- TRUE
forup[which(forup$name %in% has_SV), 
    "olSV"] <- TRUE


for (k in 2:6) {
    forup[, k] <- as.integer(forup[, k])
}

# create an upset plot from forup
cat("create upset plot\n")
pdf(sprintf("%s/upset_plot.pdf", outDir), width = 16, height = 8)
upset(
    forup, 
    sets = rev(c("Overlaps_FetalCB_Enhancer", "ABC_gene", 
             "ABC_IsNeurodev", "ABC_IsG34gene", "olSV")),
    keep.order = TRUE,         
    order.by = "freq",
    mainbar.y.label = "Number of unique\noverlapping DMRs",
    sets.x.label = "Number of annotations",
    text.scale=3.5,
    mb.ratio = c(0.7,0.3),
    point.size = 4,
    line.size = 2
)
dev.off()

cat("Get DMRs overlapping HARs\n")

harGR <- getHARs()
x <- subsetByOverlaps(dmrGR, harGR)
    cat(sprintf("DMRs overlapping HARs: %d\n", length(x)))
    

# get nearest gene, add as a column
x <- getNearestGene(x, gene_types="protein_coding")
harDMRs <- data.frame(
        chr = seqnames(x),
        start = start(x),
        end = end(x),
        nearestGene = x$nearestGene
    )

write.table(
        harDMRs, 
        file = sprintf("%s/DMRs_overlapping_HARs_%s.tsv", outDir, dt), 
        sep = "\t", 
        row.names = FALSE, 
        quote = FALSE
)


#tss <- get_gencode_anno(gene_type="protein_coding", region="tss")
}, error = function(e) {
    print(e)
}, finally={
    cat("Closing log.\n")
    sink(NULL)
})
