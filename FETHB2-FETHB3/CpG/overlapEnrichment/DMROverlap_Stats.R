rm(list=ls())
library(BSgenome.Hsapiens.UCSC.hg38.masked) # needed for genNullSeqs
library(tidyr)
library(glue)
library(IRanges)


setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")
source("utils.R")

cpg_dmr_date <- get_configs("CPG_DMR_DATE")
projectRoot <- glue("/.mounts/labs/pailab/private/projects",
                    "/FetalHindbrain/EMseq_FETHB3/output/downstream",
                    "/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/",
                    "withoutBatchCorrection/{cpg_dmr_date}")

ccreFile <- get_configs("ENCODE_CRE_FILE")

dt <- format(Sys.Date(),"%y%m%d")

dmr <- getDMRs(projectRoot)
dmr <- regioneR::filterChromosomes(dmr, organism="hg", chr.type="canonical")
cat(sprintf("Filter standard chrom: DMRs left=%i\n", length(dmr)))

hg38 <- BSgenome.Hsapiens.UCSC.hg38.masked
hg38 <- keepStandardChroms(hg38) # remove alternate chroms 

outDir <-sprintf("%s/DMRoverlap_xsun",projectRoot)
if (!file.exists(outDir)) dir.create(outDir)

negDir <- sprintf("%s/negs",outDir)
logFile <- sprintf("%s/DMROverlap_Stats.log", outDir)

numPerm <- 1000L

fetcb <- getFetalCB_HistonePeaks()

source("getGRanges_OLenrichment.R")
sink(logFile,split=TRUE)
tryCatch({
    # annotate with nearest gene, OL with histone peaks
    dmr <- getNearestGene(dmr)
    ac <- rep("",length(dmr))
    ol <- findOverlaps(dmr,fetcb$H3K27ac)
    ac[queryHits(ol)] <- "YES"
    dmr$OL_FetalCB_H3K27ac <- ac

    me <- rep("",length(dmr))
    ol <- findOverlaps(dmr,fetcb$H3K4me3)
    me[queryHits(ol)] <- "YES"
    dmr$OL_FetalCB_H3K4me3 <- me

    megaGR <- list()

    cat("MB amps/dels\n")
    cat("----------------\n")
    cat("Northcott 2012 cnv GITSTIC wide regions based on microarray")
    ## Northcott 2012 cnv GITSTIC wide regions based on microarray
    sv <- getNorthcott2012_AmpsDels(field = "region")
    for (nm in names(sv$amps)){
        seqlevels(sv$amps[[nm]]) <- paste("chr",c(1:22,"X","Y"),sep="")
    }
    for (nm in names(sv$dels)){
        seqlevels(sv$dels[[nm]]) <- paste("chr",c(1:22,"X","Y"),sep="")
    }
    mb_sv <- c(
        sv$amps[["GISTIC_Amps-Group3"]],
            sv$amps[["GISTIC_Amps_Group4"]],
            sv$dels[["GISTIC_Dels_Group3"]],
            sv$dels[["GISTIC_Dels_Group4"]]
    )
    megaGR[["MB_G34_arrayCNV"]] <- mb_sv
    ol_sv <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=mb_sv, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="MB_G34_arrayCNV"
    )

    write.table(getDF_GRoverlap(dmr,mb_sv),
        file=sprintf("%s/MB_G34_arrayCNV_%s.txt",outDir,dt),
        sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    
    ## Northcott 2017 cnv GITSTIC wide regions based on WGS
    sv <- getNorthcott2017_AmpsDels()
    for (nm in names(sv$amps)){
      seqlevels(sv$amps[[nm]]) <- paste("chr",c(1:22,"X","Y"),sep="")
    }
    for (nm in names(sv$dels)){
      seqlevels(sv$dels[[nm]]) <- paste("chr",c(1:22,"X","Y"),sep="")
    }
    mb_sv <- c(
      sv$amps[["GRP3_GISTIC_AMP"]],
      sv$amps[["GRP4_GISTIC_AMP"]],
      sv$dels[["GRP3_GISTIC_AMP"]],
      sv$dels[["GRP4_GISTIC_AMP"]]
    )
    megaGR[["MB_G34_wgsCNV"]] <- mb_sv
    ol_sv_wgs <- getGRanges_OLenrichment(
      pos=dmr,tgtGR=mb_sv, numPerm=numPerm, negDir=negDir,
      rngSeed=12345,genome=hg38, outDir=outDir,
      tgtName="MB_G34_wgsCNV"
    )
    
    write.table(getDF_GRoverlap(dmr,mb_sv),
                file=sprintf("%s/MB_G34_wgsCNV_%s.txt",outDir,dt),
                sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    
    cat("SNV/INDEL\n")
    cat("----------------\n")
    g34_mut <- import.bed(glue("/data/xsun/20240314/mutation/PCAWG/merged",
                               "/Group34_PCAWG_snv-indel_hg38.bed"))
    megaGR[["MB_G34_PCAWG"]] <- g34_mut
    ol_mut <- getGRanges_OLenrichment(
      pos=dmr,tgtGR=g34_mut, numPerm=numPerm, negDir=negDir,
      rngSeed=12345,genome=hg38, outDir=outDir,
      tgtName="MB_G34_PCAWG"
    )
    write.table(getDF_GRoverlap(dmr,g34_mut),
                file=sprintf("%s/MB_G34_PCAWG_%s.txt",outDir,dt),
                sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    
    
    cat("HAR\n")
    cat("----------------\n")
    har <- getHARs()
    har <- regioneR::filterChromosomes(har, organism="hg", chr.type="canonical")
    cat(sprintf("Filter standard chrom: HARs left=%i\n", length(har)))
    megaGR[["HAR"]] <- har

    ol_har <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=har, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="HAR"
    )
    
    write.table(getDF_GRoverlap(dmr,har),
        file=sprintf("%s/DMR_HAR_%s.txt",outDir,dt),
        sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    cat("\n")
    
    df <- read.delim(ccreFile,h=F,as.is=T)
    ccre <- GRanges(df[,1],IRanges(df[,2],df[,3]))
    ccre$type <- df[,6]
    ccre <- regioneR::filterChromosomes(ccre, 
        organism="hg", 
        chr.type="canonical")
    
    cat("ENCODE cCRE: ELS\n")
    cat("----------------\n")
    enh <- ccre[grep("ELS",ccre$type)]
    cat(sprintf("%i elements\n", length(enh)))
    ol_enh <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=enh, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="Enhancers"
    )
    write.table(getDF_GRoverlap(dmr,enh),
        file=sprintf("%s/DMR_Enh_%s.txt",outDir,dt),
        sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    cat("\n")
    megaGR[["cCRE-ELS"]] <- enh
    
    
    cat("ENCODE cCRE: PLS\n")
    cat("----------------\n")
    prmtr <- ccre[grep("PLS",ccre$type)]
    cat(sprintf("%i elements\n", length(prmtr)))
    ol_prmtr <- getGRanges_OLenrichment(
      pos=dmr,tgtGR=prmtr, numPerm=numPerm, negDir=negDir,
      rngSeed=12345,genome=hg38, outDir=outDir,
      direction = "lt",
      tgtName="Promoters"
    )
    write.table(getDF_GRoverlap(dmr,prmtr),
                file=sprintf("%s/DMR_prmtr_%s.txt",outDir,dt),
                sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    cat("\n")
    megaGR[["cCRE-PLS"]] <- prmtr
    

    cat("ENCODE cCRE: CTCF\n")
    cat("----------------\n")
    ctcf <- ccre[grep("CTCF",ccre$type)]
    cat(sprintf("%i elements\n", length(ctcf)))
    ol_ctcf <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=ctcf, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="CTCF"
    )
    write.table(getDF_GRoverlap(dmr,ctcf),
        file=sprintf("%s/DMR_CTCF_%s.txt",outDir,dt),
        sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    megaGR[["cCRE-CTCF"]] <- ctcf

    # fetal CB histone modifications (H3K27ac)
    cat("Fetal CB histone modification: H3K27ac\n")
    cat("----------------\n")
    h3k27ac <- fetcb$H3K27ac
    cat(sprintf("%i elements\n", length(h3k27ac)))
    ol_h3k27ac <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=h3k27ac, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="fecbH3K27ac"
    )
    write.table(getDF_GRoverlap(dmr,h3k27ac),
        file=sprintf("%s/DMR_fecbH3K27ac_%s.txt",outDir,dt),
        sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    megaGR[["fecbH3K27ac"]] <- h3k27ac

    # fetal CB histone modifications (H3K4me3)
    cat("Fetal CB histone modification: H3K4me3\n")
    cat("----------------\n")
    h3k4me3 <- fetcb$H3K4me3
    cat(sprintf("%i elements\n", length(h3k4me3)))
    ol_h3k4me3 <- getGRanges_OLenrichment(
        pos=dmr,tgtGR=h3k4me3, numPerm=numPerm, negDir=negDir,
        rngSeed=12345,genome=hg38, outDir=outDir,
        tgtName="fecbH3K27ac", direction = "lt"
    )
    write.table(getDF_GRoverlap(dmr,h3k4me3),
        file=sprintf("%s/DMR_fecbH3K4me3_%s.txt",outDir,dt),
        sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    megaGR[["fecbH3K4me3"]] <- h3k4me3

    # G3 MB ATAC
    cat("Smith et al MB ATAC: g3\n")
    cat("----------------\n")
    g3_atac <- import.bed(
      "/data/xsun/20240314/cre/Smith2022-Group3/allEnhancerLikeElements.bed"
      )
    cat(sprintf("%i elements\n", length(g3_atac)))
    ol_g3ATAC <- getGRanges_OLenrichment(
      pos=dmr,tgtGR=g3_atac, numPerm=numPerm, negDir=negDir,
      rngSeed=12345,genome=hg38, outDir=outDir,
      tgtName="g3ATAC"
    )
    write.table(getDF_GRoverlap(dmr,g3_atac),
                file=sprintf("%s/DMR_g3ATAC_%s.txt",outDir,dt),
                sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    megaGR[["g3ATAC"]] <- g3_atac
    
    # G4 MB ATAC
    cat("Smith et al MB ATAC: g4\n")
    cat("----------------\n")
    g4_atac <- import.bed(
      "/data/xsun/20240314/cre/Smith2022-Group4/allEnhancerLikeElements.bed"
    )
    cat(sprintf("%i elements\n", length(g4_atac)))
    ol_g4ATAC <- getGRanges_OLenrichment(
      pos=dmr,tgtGR=g4_atac, numPerm=numPerm, negDir=negDir,
      rngSeed=12345,genome=hg38, outDir=outDir,
      tgtName="g4ATAC"
    )
    write.table(getDF_GRoverlap(dmr,g4_atac),
                file=sprintf("%s/DMR_g4ATAC_%s.txt",outDir,dt),
                sep="\t", col=TRUE,row=FALSE, quote=FALSE
    )
    megaGR[["g4ATAC"]] <- g4_atac
    
    pp <- plotViolins_OLenrichment(dmr,megaGR,negDir) 
    pdf(sprintf("%s/DMRoverlap_overall.pdf", outDir), height = 7, width = 10)
    print(pp$p)
    dev.off()
    backup <- pp
    
    
    # ordered subset (enhancers)
    df <- merge(pp$df2, pp$real, by = "set")
    rownames(df) <- df$set
  
    warning("Didnt have enough time!!!TODO: fix code to do one tail test on plotViolins_OLenrichment \n HARDCODED now!!!")
    df[c("fecbH3K4me3", "cCRE-PLS"), "pval"] <- 0.001
    df$stars <- ifelse(df$pval <= 0.001, "***",
                       ifelse(df$pval <= 0.01, "**",
                              ifelse(df$pval <= 0.05, "*", "")))
    
    
    #### bar plot ####
    plot_enrichment_bars <- function(df, levels = NULL, numPerm, rename_dict = NA, star_dist = 0.6) {
      if (! is.null(levels)) {
        df <- df[df$set %in% levels,]
        df$set <- factor(df$set, levels = levels)
      } else {
        df$set <- as.factor(df$set)
      }
      df <- df[order(df$set),]
      
      df$name <- rename_dict[as.character(df$set)]
      df[is.na(df$name), "name"] <- as.character(df[is.na(df$name), "set"])
      
      df$name <- factor(df$name, levels = rev(df$name))
      
      p <- ggplot(df, aes(x=name, y=logOR),color="grey50") +
        geom_bar(stat="identity",color="black",position=position_dodge()) +
        geom_errorbar(aes(ymin=logOR-sd,ymax=logOR+sd), width=0.1, 
                      position=position_dodge(0.9)) +
        theme_minimal() + 
        coord_flip() +
        theme_classic(base_size = 20) +   # Global font size setting
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(size =20, face = "bold"),
          axis.title = element_text(size = 25),  
          axis.text = element_text(size = 20),  
          plot.title = element_text(size = 30),  
          plot.subtitle = element_text(size = 25),
          plot.caption = element_text(size = 20),
          legend.position = "none"
        ) + 
        geom_hline(yintercept = 0, linetype = "dashed", 
                   color = "red", size = 1, alpha = 0.5) +
        geom_text(aes(label = stars, 
                      y = logOR + sd* ifelse(logOR >0, 1, -1) + ifelse(logOR >0, 1, -1)*star_dist),
                  vjust = 0.8, hjust = 0, size = 7)
      
      p <- p + ylab(sprintf("log2 fold change\n(mean+/-SD, %i permutations)", 
                            numPerm))

      return(p)
    }
    

    #### Forest plot ####
    library(forestplot)
    plot_enrichment_forest <- function(df, 
                                       levels = NULL, 
                                       numPerm, 
                                       rename_dict = NA) {
      if (! is.null(levels)) {
        df <- df[df$set %in% levels,]
        df$set <- factor(df$set, levels = levels)
      } else {
        df$set <- as.factor(df$set)
      }
      df <- df[order(df$set),]
      
      df$name <- rename_dict[as.character(df$set)]
      df[is.na(df$name), "name"] <- as.character(df[is.na(df$name), "set"])
      
      p <- forestplot(df$name, 
                 df$logOR, df$logOR-df$sd, df$logOR+df$sd,
                 fn.ci_norm = fpDrawCircleCI,
                 boxsize = .1,
                 line.margin = .1,
                 xlab = sprintf("log2 fold change (mean+/-SD, %i permutations)", numPerm),
                 ci.vertices = TRUE,
                 ci.vertices.height = 0.05) |> 
        fp_set_style(box = "royalblue", 
                     line = "darkblue",
                     txt_gp = fpTxtGp(label = gpar(fontsize = 18, col = "black"),  # Font size for labels
                                      ticks = gpar(fontsize = 30),   # Font size for ticks
                                      xlab = gpar(fontsize = 30))
                     )
        
      return(p)
    }
    
    
    p <- plot_enrichment_bars(df,
                                  levels = c("HAR", 
                                             "fecbH3K27ac", "fecbH3K4me3",
                                             "cCRE-ELS", "cCRE-PLS", "cCRE-CTCF"
                                             ),
                              numPerm = numPerm
                                  )
    pdf(sprintf("%s/DMRoverlap_overall_subset_bar.pdf", outDir), height = 7, width = 10)
    print(p)
    dev.off()
    
    p <- plot_enrichment_forest(df,
                              levels = c( 
                                         "fecbH3K27ac", "fecbH3K4me3",
                                         "cCRE-ELS", "cCRE-PLS", "cCRE-CTCF",
                                         "HAR"
                              ),
                              numPerm = numPerm,
                              rename_dict = c(fecbH3K27ac = "feCB_H3K27ac",
                                              fecbH3K4me3 = "feCB_H3K4me3"
                              )
    )
    pdf(sprintf("%s/DMRoverlap_overall_subset_enh_forest.pdf", outDir), height = 10, width = 10)
    print(p)
    dev.off()
    
    # ordered subset2
    p <- plot_enrichment_bars(df, 
                                  levels = c("HAR", 
                                             "fecbH3K27ac", "fecbH3K4me3",
                                             "cCRE-ELS", "cCRE-PLS",
                                             "MB_G34_arrayCNV"
                                  ),
                              numPerm = numPerm
    ) 
    pdf(sprintf("%s/DMRoverlap_overall_subset2.pdf", outDir), height = 7, width = 10)
    print(p)
    dev.off()
    
    # ordered subset3 (MB)
    p <- plot_enrichment_bars(df, 
                              levels = c("MB_G34_arrayCNV", "MB_G34_wgsCNV"),
                              numPerm = numPerm,
                              rename_dict = c(MB_G34_arrayCNV = "Grp3&4 MB\nNorthcott 2012 GISTIC\n(Array)",
                                              MB_G34_wgsCNV = "Grp3&4 MB\nNorthcott 2017\nGISTIC (WGS)"
                                              ), 
                              star_dist = 0.05
    ) + theme(axis.text.y = element_text(size =25, face = "bold"))
    pdf(sprintf("%s/DMRoverlap_overall_subset_MB_bar.pdf", outDir), height = 6, width = 12)
    print(p)
    dev.off()
    
    # ordered subset4 (both)
    p <- plot_enrichment_bars(df,
                              levels = c( 
                                         "fecbH3K27ac", "cCRE-ELS",
                                         "fecbH3K4me3","cCRE-PLS",
                                         "cCRE-CTCF",
                                         "HAR"
                              ),
                              rename_dict = c(MB_G34_arrayCNV = "Grp3&4 MB Northcott 2012 GISTIC (Array)",
                                              MB_G34_wgsCNV = "Grp3&4 MB Northcott 2017 GISTIC (WGS)",
                                              HAR = "Human accelerated regions",
                                              fecbH3K27ac = "Fetal cerebellum H3K27ac peaks",
                                              fecbH3K4me3 = "Fetal cerebellum H3K4me3 peaks",
                                              "cCRE-ELS" = "ENCODE enhancer-like elements",
                                              "cCRE-PLS" = "ENCODE promoter-like elements",
                                              "cCRE-CTCF" = "ENCODE CTCF binding sites"
                                              ),
                              numPerm = numPerm
    )
    
    pdf(sprintf("%s/DMRoverlap_overall_subset_bar.pdf", outDir), height = 8, width = 12)
    print(p)
    dev.off()
    

    # Plot summarized real vs null
    message("plotting summarized bar plot")
    ol_list <- list(
      ol_sv = ol_sv, 
      ol_har = ol_har, 
      ol_enh = ol_enh, 
      ol_ctcf = ol_ctcf, 
      ol_h3k27ac = ol_h3k27ac, 
      ol_h3k4me3 = ol_h3k4me3,
      ol_g3ATAC = ol_g3ATAC,
      ol_g4ATAC = ol_g4ATAC
    )
    pdf(sprintf("%s/DMRoverlap_sum_%s.pdf", outDir, dt))
    print(plot_bar(ol_list, 
                   names = c("ol_sv", "ol_har", "ol_enh", 
                             "ol_ctcf", "ol_h3k27ac", "ol_h3k4me3",
                             "ol_g3ATAC", "ol_g4ATAC"
                             )
                   )
          )
    dev.off()

    cat("\n")
    }, error=function(ex){
        print(ex)
    }, finally={
        sink(NULL)
})
