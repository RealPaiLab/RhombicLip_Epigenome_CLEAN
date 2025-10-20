#' LiftOver GRanges Object using UCSC liftOver cmd tool
#' @param x (GRanges) The target intervals
#' @param chain (character) The chain file for lift over. Default hg19 to hg38.
#' @param liftOver_bin (character) The liftOver binary cmd tool directory.
#' @return (GRanges) Converted intervals. 
cmd_liftOver <- function(
    x, 
    chain = "/home/rstudio/isilon/src/ucsc-tools/chain_files/hg19ToHg38.over.chain",
    liftOver_bin = "/home/rstudio/isilon/private/xsun/Software/liftOver"
    ) {
  message(sprintf("Preparing to liftOver %d intervals.", length(x)))
  x_bed <- tempfile(fileext = ".bed")
  lifted <- tempfile(fileext = ".bed")
  unlifted <- tempfile(fileext = ".unlifted")
  cmd <- paste(liftOver_bin, x_bed, chain, lifted, unlifted)
  
  tryCatch({
    rtracklayer::export.bed(x, x_bed)
    
    message(cmd)
    system(cmd)
    
    tmp <- import.bed(lifted)
    message(sprintf("liftOver of %d intervals were successful.", length(tmp)))
  },
  error = function(e) {message(e)},
  warning = function(w) {message(w)},
  finally = {unlink(c(x_bed, lifted, unlifted))}
  ) 

  return(tmp)  
}


getGC <- function(gr){
    freqs <- alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38,gr))
    gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
    return(gc)
}

#' return BSgenome object with only autosomes, sex chromosomes, and mitochondrial.
#' Remove all other alternate chromosomes
#' @param genome (BSgenome) BSgenome object to alter
#' @return BSgenome object
keepStandardChroms <- function(genome){
    seqnames <- paste0("chr",c(1:22,"X","Y","M"))
    stopifnot(all(seqnames %in% seqnames(genome)))
    genome@user_seqnames <- setNames(seqnames, seqnames)
    genome@seqinfo <- genome@seqinfo[seqnames]
    genome
}

#' compare length, GC content and mappability of two sets of GRanges
#' @param gr1 (GRanges) interval 1
#' @param gr2 (GRanges) interval 2
#' @return (list) of ggplot objects showing width and GC content comparison
compareGR <- function(gr1, gr2){
    df1 <- data.frame(
        set=rep("set1", length(gr1)),
        wd=width(gr1),
        gc=getGC(gr1),
        chrom=seqnames(gr1)
    )
    df2 <- data.frame(
        set=rep("set2",length(gr2)),
        wd=width(gr2),
        gc=getGC(gr2),
        chrom=seqnames(gr2)
    )
    
    both <- rbind(df1,df2); 
    both$set <- as.factor(both$set)
    library(ggplot2)
    p1 <- ggplot(both, aes(x=wd,color=set))+
        geom_density() + 
        ggtitle(sprintf("length (N1=%i, N2=%i)",length(gr1),length(gr2)))
    
    p2 <- ggplot(both, aes(x=gc,color=set)) +
        geom_density()+
        ggtitle("Percent GC")

    p <- ggarrange(plotlist=list(p1=p1,p2=p2),nrow=2)
    return(p)
}


#' get RL VZ vs SVZ DMRs from Xinghan's anlaysis
getDMRs <- function(inDir){
    dmrFile <-  sprintf("%s/DMRs.csv", inDir)
      message("reading DMR")
    dmr <- read.delim(dmrFile,sep="\t",h=T,as.is=T)
    cat(sprintf("%i DMR read\n", nrow(dmr)))
    dmr <- GRanges(dmr[,1],IRanges(dmr[,2],dmr[,3]))
    return(dmr)
}


#' get HARs from Pollard lab.
getHARs <- function() {
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)

    HAR_hg19 <- "/home/rstudio/isilon/src/evolution/PollardLab_HARs/nchaes_merged_hg19.bed"
    hg19_to_hg38 <- "/home/rstudio/isilon/src/ucsc-tools/chain_files/hg19ToHg38.over.chain"

    message("reading HAR")
        har <- read.delim(HAR_hg19,sep="\t",h=F,as.is=T)
        cat(sprintf("%i HAR read\n", nrow(har)))
        har <- GRanges(har[,1],IRanges(har[,2],har[,3]),
            name=har[,4])
        har38 <- cmd_liftOver(har, hg19_to_hg38)

        cat(sprintf("%i converted to hg38\n",length(har38)))
        har38$GC <- getGC(har38)
        har38$len <- log10(width(har38))
        return(har38)
}

#' Gets amplifications/deletions from Northcott 2017 Nature (WGS based)
#' liftOver from hg19 to hg38
getNorthcott2017_AmpsDels <- function() {
  require(readxl)
  dir <- "/home/rstudio/isilon/src/MB_genomics/WGS/Northcott2017/41586_2017_BFnature22973_MOESM2_ESM.xlsx"
  
  sheetNames <- apply(expand.grid(c("GRP3", "GRP4", "SHH"), "GISTIC", c("AMP", "DEL")), 
        1, 
        paste, 
        collapse = "_")
  
  res <- list(amps = NULL, dels = NULL)
  for (sheet in sheetNames) {
    message(glue("# Processing {sheet}"))
    
    tmp <- read_excel(dir, sheet = sheet, 
                      skip = 3, n_max = 1, 
                      col_names = F, trim_ws = T
                      )
    regions <- as.vector(unlist(tmp[,2:ncol(tmp)]))
    
    message(sprintf("-- %d wide regions\n", length(regions)))
    
    regions_df <- as.data.frame(t(as.data.frame(data.frame(strsplit(regions, ":|-")))))
    colnames(regions_df) <- c("seqnames", "start", "end")
    rownames(regions_df) <- NULL
    regions_gr <- GenomicRanges::makeGRangesFromDataFrame(regions_df)
    regions_gr$name <- sprintf("%s:%d-%d",
                               seqnames(regions_gr), 
                               start(regions_gr), 
                               end(regions_gr)
                               )
    regions_gr_hg38 <- cmd_liftOver(regions_gr)
    regions_gr_hg38$source <- sheet
    
    if (grepl("AMP", sheet)) {
      res$amps[[sheet]] <- regions_gr_hg38
    } else if (grepl("DEL", sheet)) {
      res$dels[[sheet]] <- regions_gr_hg38
    } else {
      stop("Double check the sheet name")
    }
  }
  
  return(res)
}


#' Gets amplifications/deletions from Northcott 2012 Nature
#' @param field (characters) "peak" or "region" of GISTIC2 output. Default peak.
#' @details Converts hg18 coords to hg38 and returns as a list.
getNorthcott2012_AmpsDels <- function(field = "peak", verbose=FALSE){

    require(readxl)
    require(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)

    hg38 <- BSgenome.Hsapiens.UCSC.hg38

    n2012Dir <- "/home/rstudio/isilon/src/MB_genomics/SNParrays/Northcott_2012/Northcott2012_supp/nature11327-s2"
    hg18_to_hg38 <- "/home/rstudio/isilon/src/ucsc-tools/chain_files/hg18ToHg38.over.chain"

    ampFile <- sprintf("%s/2012-01-00811C-SupplementaryTable-4-GISTIC_Amps.xlsx",
        n2012Dir)
    delFile <- sprintf("%s/2012-01-00811C-SupplementaryTable-5-GISTIC_Dels.xlsx",
        n2012Dir)

    ampList <- list()
    delList <- list()

    cat("*** Amplifications in MB ***\n")
    for (sh in paste("GISTIC_Amps", c("_MB","_SHH","-Group3","_Group4"),sep="")){
        if (verbose) cat(sprintf("Reading %s\n",sh))
        if (any(grep("MB", sh))) skip <- 1 else skip <- 0
        amps <- read_excel(ampFile, sheet=sh,skip=skip)
        amps <- as.data.frame(amps)

        amps$chromosome <- paste("chr",amps$chromosome,sep="")
        
        if (field == "peak") {
          hg18 <- GRanges(amps$chromosome, IRanges(amps$peak_start, amps$peak_end))
        } else if (field == "region") {
          hg18 <- GRanges(amps$chromosome, IRanges(amps$region_start, amps$region_end))
        } else {
          stop("[field] paramter should be either peak or region.")
        }
        
        # name always use peak start and end to be consistent
        hg18$name <- paste0(amps$chromosome, ":",
                            amps$peak_start, "-",
                            amps$peak_end)
        
        if (verbose) cat(sprintf("Read %i %s(s)",length(hg18), field))
      ##  print(table(seqnames(hg18)))
        if (verbose) cat("LiftOver hg18 to hg38...")
        amps <- cmd_liftOver(hg18, hg18_to_hg38)
        amps$source <- sh
        #ch <- import.chain(hg18_to_hg38)
        #amps <- liftOver(hg18, ch)
        #if (verbose) cat(sprintf("%i ranges converted\n",length(amps)))
        #rm(hg18)
        #toomany <- which(unlist(lapply(amps,length))>1)
        #if (any(toomany)) {
        #    if (verbose) cat(sprintf("removing %i ranges with one-to-many matches\n", 
        #        length(toomany)))
        #    amps <- amps[-toomany]
        #    cat(sprintf("%i amps left\n", length(amps)))
        #}
        #amps <- unlist(amps)
    ###    amps$score <- 0
    ###    amps <- sortSeqlevels(amps)
    ###    tmp <- seqlengths(hg38);
    ###    tmp <- subset(tmp, names(tmp) %in% seqlevels(amps))
    ###    seqlengths(amps) <- tmp
        ampList[[sh]] <- amps
    }

    if (verbose) cat("\n\n*** Deletions in MB ***\n")
    for (sh in paste("GISTIC_Dels_", c("MB","SHH","Group3","Group4"),sep="")){
        if (any(grep("MB", sh))) skip <- 1 else skip <- 0
        dels <- read_excel(delFile, sheet=sh,skip=skip)
        dels <- as.data.frame(dels)

        dels$chromosome <- paste("chr",dels$chromosome,sep="")

        if (field == "peak") {
          hg18 <- GRanges(dels$chromosome, IRanges(dels$peak_start, dels$peak_end))
        } else if (field == "region") {
          hg18 <- GRanges(dels$chromosome, IRanges(dels$region_start, dels$region_end))
        } else {
          stop("[field] paramter should be either peak or region.")
        }
        
        # name always use peak start and end to be consistent
        hg18$name <- paste0(dels$chromosome, ":",
                            dels$peak_start, "-",
                            dels$peak_end)
        
        if (verbose) cat(sprintf("Read %i peaks",length(hg18)))
        ##print(table(seqnames(hg18)))
        if (verbose) cat("LiftOver hg18 to hg38...")
        dels <- cmd_liftOver(hg18, hg18_to_hg38)
        dels$source <- sh
        #ch <- import.chain(hg18_to_hg38)
        #dels <- liftOver(hg18, ch)
        #if (verbose) cat(sprintf("%i ranges converted\n",length(dels)))
        #toomany <- which(unlist(lapply(dels,length))>1)
        #if (any(toomany)) {
        #    if (verbose) cat(sprintf("removing %i ranges with one-to-many matches\n", 
        #        length(toomany)))
        #    dels <- dels[-toomany]
        #    cat(sprintf("%i dels left\n", length(dels)))
        #}
        #dels <- unlist(dels)
    ###    dels <- sortSeqlevels(dels)
    ###    dels$score <- 0
    ###    tmp <- seqlengths(hg38);
    ###    tmp <- subset(tmp, names(tmp) %in% seqlevels(dels))
    ###    seqlengths(dels) <- tmp
        delList[[sh]] <- dels
    }

    return(list(amps=ampList, dels=delList))
}


getPhastCons <- function(){
    annoRoot <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/"
    # processed in this code
    pconsFile <- sprintf("%s/phastConsElements470way.StandardChroms.withGCandLen.bed", annoRoot)

      message("reading phastCons\n")
    if (!file.exists(pconsFile)) {
        pcons <- read.delim(phastCons,h=F,as.is=T)
        pcons <- GRanges(pcons[,1],IRanges(pcons[,2],pcons[,3]))
        pcons <- keepStandardChromosomes(pcons,pruning.mode="coarse")
        cat(sprintf("%i ranges found\n", length(pcons)))
        pcons$GC <- getGC(pcons)
        pcons$len <- log10(width(pcons))
        write.table(as.data.frame(pcons), file=pconsFile,sep="\t",col=T,row=F,quote=F)
    } else {
        cat("processed file exists, reading\n")
        pconsdf <- read.delim(pconsFile,sep="\t",h=T,as.is=T)
        pcons <- GRanges(pconsdf[,1],IRanges(pconsdf[,2],pconsdf[,3]))
        pcons$GC <- pconsdf$GC
        pcons$len <- pconsdf$len
        cat(sprintf("%i ranges found\n", length(pcons)))
    }
    pcons
}

#' get nearest gene to all ranges in gr
getNearestGene <- function(gr, gene_types = NULL){
    geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v42.basic.annotation.gtf"
    genes <- rtracklayer::readGFF(geneFile)
    genes <- subset(genes, type == "gene")
    if (! is.null(gene_types)) {
      genes <- subset(genes, gene_type %in% gene_types)
    }

    genes$TSS <- genes$start
    genes$TSS[which(genes$strand=="-")] <- genes$end[which(genes$strand=="-")]

    geneGR <- GRanges(genes$seqid, 
        IRanges(genes$start, genes$end),
        name=genes$gene_name
    ) 

    n <- nearest(gr, geneGR)
    gr$nearestGene <- geneGR$name[n] 

    gr
}

#' get fetal cerebellum histone peaks as GRanges
getFetalCB_HistonePeaks <- function() {
    pDir <- "/home/rstudio/isilon/private/xsun/output/ncMutMB/20240314/cre/Aldinger-FetalCB/raw"

    ac1 <- import.bed(sprintf("%s/1_27907-102M_H3K27Ac_hg38.bed",pDir))
    ac2 <- import.bed(sprintf("%s/2_27556-132M_H3K27Ac_hg38.bed",pDir))

    me1 <- import.bed(sprintf("%s/3_27907-102M_H3K4me3_hg38.bed",pDir))
    me2 <- import.bed(sprintf("%s/4_27556-132M_H3K4me3_hg38.bed",pDir))

    return(list(H3K27ac=IRanges::intersect(ac1,ac2), 
                H3K4me3=IRanges::intersect(me1,me2),
                H3K27ac_union=IRanges::union(ac1,ac2)
                ))
}



#' plot bar plot to summarize the real and null permutation results
#' @param ol_list (list) a list of overlap output from plotGRanges_OLenrichment function
#' @param names (vector) a vector of names of each overlap
#' @return (ggplot) bar plot
plot_bar <- function(ol_list, names) {
  summarize_ol <- function(ol) {
    return(data.frame( 
                      median = c(ol$overlap_pos, median(ol$overlap_negs)),
                      min = c(NA, min(ol$overlap_negs)), 
                      max = c(NA, max(ol$overlap_negs)),
                      cat = c("real", "null"), 
                      pval = rep(ol$pval, 2)
                      ))
  } 
  
  summarized_df <- Reduce(rbind, lapply(ol_list, summarize_ol))
  summarized_df$name <- rep(names, each = 2)
  
  p_labels <- na.omit(summarized_df) %>%
    mutate(label = paste0("p = ", formatC(pval, format = "e", digits = 2)))
  
  
  p <- ggplot(summarized_df, aes(x = cat, y = median)) + 
    geom_bar(stat = "identity", fill = "skyblue", alpha = 0.8, position = "dodge") + 
    facet_wrap(~ name) + 
    geom_errorbar(aes(x=cat, ymin=min, ymax=max), width=0.4, colour="orange", alpha=0.9) +
    theme_minimal()
  
  p <- p + geom_text(data = p_labels, aes(x = 1.5, y = max(summarized_df$median, na.rm = T), label = label),
                vjust = -0.5, size = 3, color = "red", inherit.aes = FALSE)
  
  return(p)
}
