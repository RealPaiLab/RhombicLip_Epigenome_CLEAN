# Collect the Cicero results from the Mannens project and analyse
rm(list=ls())

library(Seurat)
library(ggplot2)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(Signac)

inDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024_Cicero/250611"
geneDef <- "/home/rstudio/isilon/src/gencode/GRCh38/gencode.v42.basic.annotation.gtf"
inFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024/250610/Mannens_2024_seurat.qs"

outDir <- inDir

cat("Compiling all Cicero results...\n")
fList <- list.files(inDir, pattern = "cicero", full.names = TRUE)
# read in cicero results as a list and combine
t0 <- Sys.time()
ciceroList <- lapply(fList, function(x) {
    print(basename(x))
  load(x)
  conns
})
print(Sys.time()-t0)
# combine the list into a single data frame
cat("Combining Cicero results into a dataframe\n")
t0 <- Sys.time()
ciceroDF <- do.call(rbind, ciceroList)
print(Sys.time()-t0)

cat("Writing combined results to file\n")
outFile <- sprintf("%s/Mannens_cicero_conns.txt", outDir)
write.table(ciceroDF, file=outFile, sep="\t", 
    quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("Reading gene definition to build gene activity matrix")
gene_anno <- rtracklayer::import(geneDef)
gene_df <- as.data.frame(gene_anno)
gene_df$chromosome <- as.character(gene_df$seqnames)
gene_df$gene <- gene_df$gene_id
gene_df$transcript <- gene_df$transcript_id
gene_df$symbol <- gene_df$gene_name
cat(sprintf("Read %i genes from GTF file\n", 
    nrow(gene_df)))                                                    

cat("* Computing gene activity scores\n")
gene_anno <- gene_df
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)
gene_annotation_sub <- gene_annotation_sub[,
  c("chromosome", "start", "end", "symbol")]
names(gene_annotation_sub)[4] <- "gene"

cat("* Now reading in the Seurat object and converting to CellDataSet\n")
srat <- qs::qread(inFile)
DefaultAssay(srat) <- "peaks"
srat.cds <- as.cell_data_set(x=srat)
srat.cds <- annotate_cds_by_site(srat.cds,gene_annotation_sub)

cat("* Build gene activity matrix\n")
conns <- ciceroDF
t0 <- Sys.time()
unnorm_ga <- build_gene_activity_matrix(srat.cds,conns)
print(Sys.time()-t0)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]

# Monocle3: creates the num_genes_expressed column
srat.cds <- detect_genes(srat.cds) 

num_genes <- pData(srat.cds)$num_genes_expressed
names(num_genes) <- row.names(pData(srat.cds))

cicero_gene_activities <- normalize_gene_activities(
    unnorm_ga, num_genes)

qs::qsave(cicero_gene_activities, 
    file=sprintf("%s/Mannens_cicero_gene_activities_allchroms.qs", outDir))
qs::qsave(srat.cds, 
    file=sprintf("%s/Mannens_cicero_srat_cds_annotated_allchroms.qs", outDir))
qs::qsave(conns, 
    file=sprintf("%s/Mannens_cicero_conns_allchroms.qs", outDir))