# run cicero to compute co-accessibility scores
rm(list=ls())

library(Seurat)
library(ggplot2)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(Signac)

inFile <- "/home/rstudio/isilon/src/neurodev-genomics/multiome/Mannens_2024/from_anders/rl_micro_mannens.qs"
chromSizes <- "/home/rstudio/isilon/src/ucsc-goldenpath/hg38/hg38-20221124T1111.chrom.sizes"
geneDef <- "/home/rstudio/isilon/src/gencode/GRCh38/gencode.v42.basic.annotation.gtf"


outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Mannens2024_Cicero"
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}


logFile <- sprintf("%s/Mannens_Cicero.log", outDir)
sink(logFile,split=TRUE)
tryCatch({
t0 <- Sys.time()
srat <- qs::qread(inFile)
t1 <- Sys.time()
message(sprintf("Time taken to read file: %s", t1 - t0))

# process the ATAC peaks
DefaultAssay(srat) <- "peaks"
srat <- FindTopFeatures(srat, min.cutoff = 10)
t0 <- Sys.time()
srat <- RunTFIDF(srat)
t1 <- Sys.time()
message(sprintf("Time taken for TF-IDF: %s", t1 - t0))
t0 <- Sys.time()
srat <- RunSVD(srat)
t1 <- Sys.time()
message(sprintf("Time taken for SVD: %s", t1 - t0))
srat <- RunUMAP(srat, reduction = "lsi", dims=2:30)


DefaultAssay(srat) <- "peaks"
# running Cicero
# instructions from tutorial: https://stuartlab.org/signac/articles/cicero
srat.cds <- as.cell_data_set(x=srat)
t0 <- Sys.time()
srat.cicero <- make_cicero_cds(srat.cds, 
    reduced_coordinates=reducedDims(srat.cds)$UMAP
)
t1 <- Sys.time()
print(sprintf("Time taken to make Cicero CDS: %1.2fs",t1- t0))

chrom <- read.delim(chromSizes,header=FALSE)
#chrom <- subset(chrom, V1 %in% c("chr2","chr3","chr5"))

cat("Running Cicero for each chromosome\n")
for (i in 1:nrow(chrom)){
    chrom2 <- subset(chrom, V1 == chrom$V1[i])
    print(chrom2)
    t0 <- Sys.time()
    conns <- run_cicero(
        cds = srat.cicero,
        genomic_coords = chrom2,
        sample_num = 100
    )
    print(Sys.time()-t0)

    save(conns, file = sprintf("%s/cicero_conns_%s.RData", outDir, chrom2$V1))
    message(sprintf("Cicero connections for %s saved.", chrom2$V1))
}
cat("done running for all chromosomes\n")
browser()

conns <- list()
for (i in 1:nrow(chrom)){
    chrom2 <- subset(chrom, V1 == chrom$V1[i])
    print(chrom2)
    load(sprintf("%s/cicero_conns_%s.RData", outDir, chrom2$V1))
    conns[[as.character(chrom2$V1)]] <- conns
}
conns <- do.call(rbind, conns)
cat(sprintf("Compiled %i connections from Cicero\n", nrow(conns)))

###gene_df <- as.data.frame(gene_anno)
###gene_anno <- rtracklayer::import(geneDef)
###gene_df$chromosome <- as.character(gene_df$seqnames)
###gene_df$gene <- gene_df$gene_id
###gene_df$transcript <- gene_df$transcript_id
###gene_df$symbol <- gene_df$gene_name
####gene_df <- subset(gene_df, gene_type == "protein_coding")
###cat(sprintf("Read %i genes from GTF file\n", nrow(gene_df)))                                                    
###cat("Computing gene activity scores\n")
###cat("* Creating DF with TSS positions\n")
###
###gene_anno <- gene_df
###pos <- subset(gene_anno, strand == "+")
###pos <- pos[order(pos$start),] 
#### remove all but the first exons per transcript
###pos <- pos[!duplicated(pos$transcript),] 
#### make a 1 base pair marker of the TSS
###pos$end <- pos$start + 1 
###neg <- subset(gene_anno, strand == "-")
###neg <- neg[order(neg$start, decreasing = TRUE),] 
###
#### remove all but the first exons per transcript
###neg <- neg[!duplicated(neg$transcript),] 
###neg$start <- neg$end - 1
###gene_annotation_sub <- rbind(pos, neg)
###gene_annotation_sub <- gene_annotation_sub[,
###
###  c("chromosome", "start", "end", "symbol")]
###names(gene_annotation_sub)[4] <- "gene"
###cat("* Annotating CDS by site\n")
###t0 <- Sys.time()
###
###srat.cds <- annotate_cds_by_site(srat.cds,gene_annotation_sub)
###t1 <- Sys.time()
###message(sprintf("Time taken to annotate CDS by site: %s", t1 - t0))
###cat("* Build gene activity matrix\n")
###t0 <- Sys.time()
###
###unnorm_ga <- build_gene_activity_matrix(srat.cds,conns)
###print(Sys.time()-t0)
#### remove any rows/columns with all zeroes
###
###unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
###                       !Matrix::colSums(unnorm_ga) == 0]
###
#### make a list of num_genes_expressed
###srat.cds <- detect_genes(srat.cds) # creates the num_genes_expressed column
###num_genes <- pData(srat.cds)$num_genes_expressed
###names(num_genes) <- row.names(pData(srat.cds))
###
###cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)



}, error = function(ex) {
  print(ex)
}, finally= {
    cat("Closing log.\n")
  sink(NULL)
}
)