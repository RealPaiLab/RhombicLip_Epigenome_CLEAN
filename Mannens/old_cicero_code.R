# probably outdated and can be removed.
#### running Cicero
#### instructions from tutorial: https://stuartlab.org/signac/articles/cicero
###srat.cds <- as.cell_data_set(x=srat)
###t0 <- Sys.time()
###srat.cicero <- make_cicero_cds(srat.cds, 
###    reduced_coordinates=reducedDims(srat.cds)$UMAP
###)
###t1 <- Sys.time()
###print(sprintf("Time taken to make Cicero CDS: %1.2fs",t1- t0))
###
###chrom <- read.delim(chromSizes,header=FALSE)
####chrom <- subset(chrom, V1 %in% c("chr2","chr3","chr5"))
###
###for (i in 1:nrow(chrom)){
###conns <- run_cicero(
###    cds = srat.cicero,
###    genomic_coords = chrom2,
###    sample_num = 100
###)
###
###
#### plot distribution of coaccess as a violin and boxplot
###conns$chrom <- "chr17"
###p <- ggplot(conns, aes(x=chrom, y = coaccess)) +
###    geom_violin(fill = "lightblue", alpha = 0.5) +
###    geom_boxplot(width = 0.1, fill = "red", outlier.shape = NA) +
###    labs(title = "Distribution of Cicero Coaccess Weights",
###         x = "Coaccess Weight",
###         y = "Density") +
###    ylim(c(0,0.5)) + 
###    theme_minimal(base_size = 20)
###
###ggsave(
###    filename = sprintf("%s/cicero_weights_distribution.pdf", outDir),
###    plot = p,
###    width = 8, height = 6
###)
###
###gene_anno <- rtracklayer::import(geneDef)
###gene_df <- as.data.frame(gene_anno)
###gene_df$chromosome <- as.character(gene_df$seqnames)
###gene_df$gene <- gene_df$gene_id
###gene_df$transcript <- gene_df$transcript_id
###gene_df$symbol <- gene_df$gene_name
####gene_df <- subset(gene_df, gene_type == "protein_coding")
###cat(sprintf("Read %i genes from GTF file\n", nrow(gene_df)))                                                    
###
###cat("Computing gene activity scores\n")
###cat("* Creating DF with TSS positions\n")
###gene_anno <- gene_df
###pos <- subset(gene_anno, strand == "+")
###pos <- pos[order(pos$start),] 
#### remove all but the first exons per transcript
###pos <- pos[!duplicated(pos$transcript),] 
#### make a 1 base pair marker of the TSS
###pos$end <- pos$start + 1 
###
###neg <- subset(gene_anno, strand == "-")
###neg <- neg[order(neg$start, decreasing = TRUE),] 
#### remove all but the first exons per transcript
###neg <- neg[!duplicated(neg$transcript),] 
###neg$start <- neg$end - 1
###
###gene_annotation_sub <- rbind(pos, neg)
###gene_annotation_sub <- gene_annotation_sub[,
###  c("chromosome", "start", "end", "symbol")]
###names(gene_annotation_sub)[4] <- "gene"
###
###cat("* Annotating CDS by site\n")
###t0 <- Sys.time()
###srat.cds <- annotate_cds_by_site(srat.cds,gene_annotation_sub)
###t1 <- Sys.time()
###message(sprintf("Time taken to annotate CDS by site: %s", t1 - t0))
###
###cat("* Build gene activity matrix\n")
###t0 <- Sys.time()
###unnorm_ga <- build_gene_activity_matrix(srat.cds,conns)
###print(Sys.time()-t0)
###
#### remove any rows/columns with all zeroes
###unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
###                       !Matrix::colSums(unnorm_ga) == 0]
###
#### make a list of num_genes_expressed
###srat.cds <- detect_genes(srat.cds) # creates the num_genes_expressed column
###num_genes <- pData(srat.cds)$num_genes_expressed
###names(num_genes) <- row.names(pData(srat.cds))
###
###cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
###
#### plot view around BRCA1 gene.
###pdf(sprintf("%s/cicero_gene_activity_BRCA1.pdf", outDir), 
###  width = 10, height = 6)
###plot_connections(conns, 
###                 alpha_by_coaccess = TRUE, 
###                 "chr17", 43016427, 43153232,
###                 gene_model = gene_anno, 
###                 coaccess_cutoff = 0.2, 
###                 connection_width = .5, 
###                 collapseTranscripts = "longest" )
###dev.off()              