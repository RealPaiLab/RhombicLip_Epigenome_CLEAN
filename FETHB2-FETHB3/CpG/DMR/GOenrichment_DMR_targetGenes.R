# Enrichment of DMR target genes in rhombic lip DEGs and gene sets.

rm(list=ls())
library(gprofiler2)

abcFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/Nasser-Neuronal-ABC_creTarget_hg38.bed"

# RL-VZ vs. RL-SVZ gene expression data from Haldipur et al. 2019
phenoFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/input/Haldipur_RL_VZSVZ_counts/aldinger_rnaseq_0218_all.star_fc.metadata_rlvzsvz.txt"
rnaFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/AldingerMillen_LCM_RNAseq/input/Haldipur_RL_VZSVZ_counts/aldinger_rnaseq_0218_rl_0518.star_fc.counts_rlvzsvz.txt"

dmr2genes <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_link2Genes_ABC/251007/DMR_AnnotatedAll_251007.tsv"

outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_TargetGene_Enrichment"

outDir <- sprintf("%s/GOenrichment_%s", outDir, format(Sys.Date(),"%y%m%d"))
if (!file.exists(outDir)) dir.create(outDir, recursive = FALSE)

logFile <- sprintf("%s/GOenrichment_DMR_targetGenes_logfile.txt", outDir)
sink(logFile, append=FALSE, split=TRUE)

tryCatch({
    dmr2genes_data <- read.delim(dmr2genes, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    fg_genes <- unique(dmr2genes_data$ABC_gene)

    abc <- read.delim(abcFile, header=FALSE, stringsAsFactors = FALSE)
    abc <-  subset(abc, V1 %in% paste("chr", c(1:22, "X", "Y"),sep=""))
    bg_genes <- unique(abc[,4])

    cat(sprintf("%i FG genes and %i BG genes loaded.\n\n", 
        length(fg_genes), length(bg_genes)))

    source("../../gprofiler2_helpers.R")
        # run pathway enrichment analysis (from `gprofiler2_helpers.R`)

    gost_res <- run_gost(
          query = fg_genes,
          organism = "hsapiens",
    #      organism =  "gp__OUFl_gpIE_RyE",
          significant = TRUE,
          sources = c("GO:BP","GO:CC","REAC","WP"),
          evcodes = TRUE,
          correction_method = "fdr",
          custom_bg = bg_genes,
          filename = file.path(outDir, "enriched_pathways")
        )

    if (!is.null(gost_res) && nrow(gost_res$result) > 0) {
          res <- gost_res$result
          cat(sprintf("Found %i significantly enriched pathways.\n", nrow(res)))

          print(head(res$term_id))
          # save results table
          df <- gost_res$result
          df$parents <- NULL
          write.table(
            df,
            file = sprintf("%s/gost_enrichment_results.tsv", outDir),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
          )


        # plot a barplot of all enriched terms. length of barplot should be
        # -log10(p_value), and label should be term_name
        # use font Helvetica, size 18 for all text
        p <- ggplot(res, aes(x=reorder(term_name, -p_value), y=-log10(p_value))) +
          geom_bar(stat="identity", fill="lightblue") +
          coord_flip() +
          xlab("Enriched Terms") +
          ylab("-log10(p-value)") +
          ggtitle("Enrichment of DMR Target Genes") +
          theme_bw() 
        p <- p + theme(text=element_text(family="Helvetica", size=18))
        ggsave(filename = sprintf("%s/enriched_terms_barplot.pdf", outDir),
          plot = p, width = 8, height = 6)

        } else {
          cat("No significant enrichment found.\n")
        }
 

    # ----------------------------------------------------------------------------

    cat("Finding DEG in >14 pcw RL-VZ vs. RL-SVZ...\n")
    pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
    dat <- read.delim(rnaFile,sep="\t",h=T,as.is=T)
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

    #message("matching pheno to data")
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

    group <- factor(pheno$tissue)
    y <- DGEList(counts=dat, group=group)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y,method="TMM")
    design <- model.matrix(~group)
    y <- estimateDisp(y,design)

    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    tt <- lrt$table
    tt$FDR <- p.adjust(tt$PValue,method="BH")
    cat(sprintf("%i genes survive FDR correction\n\n", sum(tt$FDR < 0.05)))
    deg_full <- tt

    # ----------------------------------------------------------------------------
    cat("Now computing enrichment of DMR target genes in DEGs\n\n")
    common <- intersect(rownames(deg_full), bg_genes)
    cat(sprintf("Found %i common genes between DEG results and BG genes.\n", length(common)))
    deg_full <- deg_full[common, ]
    cat(sprintf("After subsetting, DEG results has %i genes.\n", nrow(deg_full)))

    cat(sprintf("Before intersecting common, found %i FG genes and %i BG genes.\n", 
        length(fg_genes), length(bg_genes)))
    fg_genes <- intersect(fg_genes, common)
    bg_genes <- intersect(bg_genes, common)
    cat(sprintf("After intersecting common, found %i FG genes and %i BG genes.\n", 
        length(fg_genes), length(bg_genes)))

    deg <- rownames(deg_full[which(deg_full$FDR < 0.05),])
    not_deg <- setdiff(rownames(deg_full), deg)
    cat(sprintf("Found %i DEGs (Q<0.05) between RL-VZ and RL-SVZ\n", length(deg)))
    cat(sprintf("Found %i non-DEGs\n", length(not_deg)))

    # Fisher's Exact Tet of FG genes being DEGs vs. BG genes being DEGs
    #            in_DEG not_in_DEG
    # FG_genes     a        b
    # BG_genes     c        d
    a <- length(intersect(fg_genes, deg))
    b <- length(intersect(fg_genes, not_deg))
    c <- length(intersect(bg_genes, deg))
    d <- length(intersect(bg_genes, not_deg))
    fisher_matrix <- matrix(c(a,b,c,d), nrow=2,
        dimnames=list(c("FG_genes","BG_genes"), c("in_DEG","not_in_DEG")))
    fisher_result <- fisher.test(fisher_matrix, alternative="greater")
    cat("Fisher's exact test results for FG genes being DEGs vs. BG genes being DEGs:\n")
    #print(fisher_result)

    cat(sprintf("%1.2f%% DMR-target genes are DEGs (n=%i out of %i)\n", 
        100* a / length(fg_genes), a, length(fg_genes)))
    cat(sprintf("%1.2f%% BG genes are DEGs (n=%i out of %i)\n\n", 
        100* c / length(bg_genes), c, length(bg_genes)))
    cat(sprintf("FET p < %1.2e\n", fisher_result$p.value))

    # plot these fractions as stacked barplot
    df_plot <- data.frame(
        GeneSet = c("DMR Target Genes", "Background Genes"),
        DEG = c(a, c),
        Not_DEG = c(b, d)
    )
    df_plot_melt <- reshape2::melt(df_plot, id.vars = "GeneSet",
        variable.name = "Status", value.name = "Count")
    p <- ggplot(df_plot_melt, aes(x=GeneSet, y=Count, fill=Status)) +
        geom_bar(stat="identity",position="fill") + 
        ggtitle("Enrichment of DMR Target Genes in RL-VZ vs. RL-SVZ DEGs") +
        ylab("Number of Genes") +
        theme_bw()
    p <- p + theme(text=element_text(family="Helvetica", size=18))
    ggsave(filename = sprintf("%s/DMR_targetGenes_in_RLVZvs
RLVZ_DEGs_barplot.pdf", outDir),
        plot = p, width = 8, height = 6)    
    

}, error=function(ex){
    print(ex)
}, finally={
    sink()
})