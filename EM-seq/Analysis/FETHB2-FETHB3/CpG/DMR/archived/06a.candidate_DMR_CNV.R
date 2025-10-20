rm(list=ls())

library(glue)

setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")
source("../../../../../MPRA/oligo_design/utils.R")
source("../overlapEnrichment/utils.R")

dt <- format(Sys.Date(),"%y%m%d")
cpg_dmr_date <- get_configs("CPG_DMR_DATE")


### load data ###
## hars ##
hars <- getHARs()


## vista ##
all_lines <- readLines(get_configs("VISTA_HB"))
filtered_lines <- grep("^>", all_lines, value = TRUE)
loc <- stringr::str_split(filtered_lines, "\\|", simplify = T)[,2]
loc <- stringr::str_trim(loc)
element <- stringr::str_split(filtered_lines, "\\|", simplify = T)[,3]
element <- sub(" ", "_", stringr::str_trim(element))
element <- glue("VISTA_{element}")

vista_hb <- as.data.frame(stringr::str_split(loc, ":|-", simplify = T))
colnames(vista_hb) <- c("seqnames", "start", "end")
vista_hb_gr_hg19 <- GenomicRanges::makeGRangesFromDataFrame(vista_hb)
vista_hb_gr_hg38 <- liftOver_gr(vista_hb_gr_hg19)
vista_hb_gr_hg38$name <- element

## tss ##
tss <- get_gencode_anno(gene_type = "protein_coding",region = "tss")


## all hypo segs ##
seg_dir <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMVs/CTsnv_excluded/240712/Segmentation"

seg_list <- list()
for (file in list.files(seg_dir, pattern = "*.bed")) {
  path <- glue("{seg_dir}/{file}")
  tmp <- import.bed(path)
  tmp <- tmp[tmp$name %in% c("2")]
  
  seg_list[[file]] <- tmp
}
narrowPeakList <- GRangesList(seg_list)
peakCov <- IRanges::coverage(narrowPeakList)
covered_ranges <- IRanges::slice(peakCov, lower=8, rangesOnly=T)
res <- GRanges(covered_ranges)



## dmr ##
#class <- "enhancer"
class <- "dmr"
if (class == "enhancer") {
  dmr <- c(res, get_cpg_dmrs())
} else if (class == "dmr") {
  dmr <- get_cpg_dmrs() 
}

enh <- c(res, get_cpg_dmrs())


## MB atac ##
g3_atac <- import.bed(
  "/data/xsun/20240314/cre/Smith2022-Group3/allEnhancerLikeElements.bed"
)
g4_atac <- import.bed(
  "/data/xsun/20240314/cre/Smith2022-Group4/allEnhancerLikeElements.bed"
)

## g34 cnv ##
all_g34_genes <- F

if (all_g34_genes) {
  g34_genes <- get_g34genes()
  g34_genes_2 <- g34_genes
} else {
  g34_genes <- c("PRDM6","RUNX1T1","CBFA2T3","KDM2B","GFI1B","GFI1","OTX2","MYC","MYCN","CDK6","FANCA", "GSE1") # only those amp/del studied, need to update to be more precise
  g34_genes_2 <- get_g34genes()
}


#g34_genes_2 <- get_g34genes()
g34_genes_2 <- c(get_g34genes(), get_neurodevGenes())

#g34_genes_2 <- unique(tss$name)
g34_tss <- tss[tss$name %in% g34_genes]
g34_tss_2 <- tss[tss$name %in% g34_genes_2]

# cnv <- getNorthcott2012_AmpsDels(field = "region")
# cnv <- getNorthcott2012_AmpsDels(field = "peak")
# 
cnv <- getNorthcott2017_AmpsDels()
g34_mut <- import.bed(glue("/data/xsun/20240314/mutation/PCAWG/merged/Group34_PCAWG_snv-indel_hg38.bed"))

#g34_cnv <- c(cnv$amps$`GISTIC_Amps-Group3`, cnv$amps$`GISTIC_Amps-Group4`,
#             cnv$dels$GISTIC_Dels_Group3, cnv$dels$GISTIC_Dels_Group4)
#amp_cnv <- c(cnv$amps$`GISTIC_Amps-Group3`, cnv$amps$GISTIC_Amps_Group4)
#amp_cnv <- subsetByOverlaps(amp_cnv, g34_tss, invert = T)
#del_cnv <- c(cnv$dels$GISTIC_Dels_Group3, cnv$dels$GISTIC_Dels_Group4)
#del_cnv <- subsetByOverlaps(del_cnv, g34_tss, invert = T)

#amp_cnv_dmr <- subsetByOverlaps(amp_cnv, IRanges::intersect(dmr, c(g3_atac, g4_atac)))

amp_cnv <- c(cnv$amps$GRP3_GISTIC_AMP, cnv$amps$GRP4_GISTIC_AMP)
del_cnv <- c(cnv$dels$GRP3_GISTIC_DEL, cnv$dels$GRP4_GISTIC_DEL)
amp_cnv_dmr <- subsetByOverlaps(amp_cnv, dmr)
del_cnv_dmr <- subsetByOverlaps(del_cnv, dmr)


cnv_wgs <- getNorthcott2017_AmpsDels()
cnv_aff <- getNorthcott2012_AmpsDels(field = "region")
amp_cnv <- c(subsetByOverlaps(cnv_wgs$amps$GRP3_GISTIC_AMP, cnv_aff$amps$`GISTIC_Amps-Group3`), 
  subsetByOverlaps(cnv_wgs$amps$GRP4_GISTIC_AMP, cnv_aff$amps$GISTIC_Amps_Group4)
  )
del_cnv <- c(subsetByOverlaps(cnv_wgs$dels$GRP3_GISTIC_DEL, cnv_aff$dels$GISTIC_Dels_Group3), 
             subsetByOverlaps(cnv_wgs$dels$GRP4_GISTIC_DEL, cnv_aff$dels$GISTIC_Dels_Group4)
)
amp_cnv_dmr <- subsetByOverlaps(amp_cnv, dmr)
del_cnv_dmr <- subsetByOverlaps(del_cnv, dmr)


#radius <- 2e6
#radius <- 5e5
#radius <- 5e6
radius <- 2e5
#radius = 1

amp_cnv_dmr_resized <- resize(amp_cnv_dmr, 
                              width = width(amp_cnv_dmr) + radius*2, 
                              fix = "center"
                              )
start(amp_cnv_dmr_resized) <- ifelse(start(amp_cnv_dmr_resized) < 0, 1, start(amp_cnv_dmr_resized))
del_cnv_dmr_resized <- resize(del_cnv_dmr, 
                              width = width(del_cnv_dmr) + radius*2, 
                              fix = "center"
)

start(del_cnv_dmr_resized) <- ifelse(start(del_cnv_dmr_resized) < 0, 1, start(del_cnv_dmr_resized))


amp_tmp <- amp_cnv_dmr
#amp_tmp$name <- NULL
amp_tmp$type <- "AMP"
amp_tmp$resized_cnv_dmr <- paste0(seqnames(amp_cnv_dmr_resized), ":",
                                  start(amp_cnv_dmr_resized), "-",
                                  end(amp_cnv_dmr_resized)
)

amp_tmp$n_ol_dmrs <- unlist(lapply(split(amp_tmp), function(x) {
  return(length(subsetByOverlaps(dmr, x)))
}))

amp_tmp$ol_genes <- unlist(lapply(split(amp_tmp), function(x) {
  return(paste(intersect(g34_genes_2, unique(subsetByOverlaps(tss, x)$name)), collapse = ", "))
}))

amp_tmp$nearby_g34_genes <- unlist(lapply(split(amp_cnv_dmr_resized), function(x) {
  return(paste(intersect(g34_genes_2, unique(subsetByOverlaps(tss, x)$name)), collapse = ", "))
}))

amp_tmp$ol_hars <- unlist(lapply(split(amp_tmp), function(x) {
  return(paste(subsetByOverlaps(hars, x)$name, collapse = ", "))
}))

amp_tmp$resized_ol_hars <- unlist(lapply(split(amp_cnv_dmr_resized), function(x) {
  return(paste(subsetByOverlaps(hars, x)$name, collapse = ", "))
}))

amp_tmp$ol_vista_hb <- unlist(lapply(split(amp_tmp), function(x) {
  return(paste(subsetByOverlaps(vista_hb_gr_hg38, x)$name, collapse = ", "))
}))

amp_tmp$resized_ol_hb <- unlist(lapply(split(amp_cnv_dmr_resized), function(x) {
  return(paste(subsetByOverlaps(vista_hb_gr_hg38, x)$name, collapse = ", "))
}))

amp_tmp$n_mut_samp <- unlist(lapply(split(amp_tmp), function(x) {
  length(unique(subsetByOverlaps(g34_mut, subsetByOverlaps(enh, x))$name))
}))

amp_tmp$dmr_len <- unlist(lapply(split(amp_tmp), function(x) {
  sum(width(unique(subsetByOverlaps(enh, x))))
}))


del_tmp <- del_cnv_dmr
del_tmp$type <- "DEL"
del_tmp$resized_cnv_dmr <- paste0(seqnames(del_cnv_dmr_resized), ":",
                                      start(del_cnv_dmr_resized), "-",
                                      end(del_cnv_dmr_resized)
)

del_tmp$n_ol_dmrs <- unlist(lapply(split(del_tmp), function(x) {
  return(length(subsetByOverlaps(dmr, x)))
}))


del_tmp$ol_genes <- unlist(lapply(split(del_tmp), function(x) {
  return(paste(intersect(g34_genes_2, unique(subsetByOverlaps(tss, x)$name)), collapse = ", "))
}))

del_tmp$nearby_g34_genes <- unlist(lapply(split(del_cnv_dmr_resized), function(x) {
  return(paste(intersect(g34_genes_2, unique(subsetByOverlaps(tss, x)$name)), collapse = ", "))
}))

del_tmp$ol_hars <- unlist(lapply(split(del_tmp), function(x) {
  return(paste(subsetByOverlaps(hars, x)$name, collapse = ", "))
}))

del_tmp$resized_ol_hars <- unlist(lapply(split(del_cnv_dmr_resized), function(x) {
  return(paste(subsetByOverlaps(hars, x)$name, collapse = ", "))
}))

del_tmp$ol_vista_hb <- unlist(lapply(split(del_tmp), function(x) {
  return(paste(subsetByOverlaps(vista_hb_gr_hg38, x)$name, collapse = ", "))
}))

del_tmp$resized_ol_hb <- unlist(lapply(split(del_cnv_dmr_resized), function(x) {
  return(paste(subsetByOverlaps(vista_hb_gr_hg38, x)$name, collapse = ", "))
}))

del_tmp$n_mut_samp <- unlist(lapply(split(del_tmp), function(x) {
  length(unique(subsetByOverlaps(g34_mut, subsetByOverlaps(enh, x))$name))
}))

del_tmp$dmr_len <- unlist(lapply(split(del_tmp), function(x) {
  sum(width(unique(subsetByOverlaps(enh, x))))
}))


combined_tmp <- c(amp_tmp, del_tmp)

combined_tmp$radius <- radius
combined_tmp$class <- class
combined_tmp$rate <- combined_tmp$n_mut_samp/combined_tmp$dmr_len

# cna count
pcawg_cna_ind <- import.bed("/data/xsun/20241005_PCAWG_MB_muts/combined_cna_hg38.bed")

combined_tmp$pcawg_dup_dmr <- unlist(lapply(split(combined_tmp), function(x) {
  dmr_ol_gistic <- subsetByOverlaps(dmr, x)
  
  c <- unlist(lapply(split(dmr_ol_gistic), function(y) {
    return(length(subsetByOverlaps(pcawg_cna_ind[pcawg_cna_ind$name>0], y)))
  })) 
  names(c) <- paste0(seqnames(dmr_ol_gistic), ":",
                     start(dmr_ol_gistic), "-",
                     end(dmr_ol_gistic)
  )
  c <- sort(c, decreasing = T)
  
  return(paste0(names(c), " (", c, ")", collapse = "; "))
}))

combined_tmp$pcawg_del_dmr <- unlist(lapply(split(combined_tmp), function(x) {
  dmr_ol_gistic <- subsetByOverlaps(dmr, x)
  
  c <- unlist(lapply(split(dmr_ol_gistic), function(y) {
    return(length(subsetByOverlaps(pcawg_cna_ind[pcawg_cna_ind$name<0], y)))
  })) 
  names(c) <- paste0(seqnames(dmr_ol_gistic), ":",
                     start(dmr_ol_gistic), "-",
                     end(dmr_ol_gistic)
  )
  c <- sort(c, decreasing = T)
  
  return(paste0(names(c), " (", c, ")", collapse = "; "))
}))


combined <- as.data.frame(combined_tmp)

#write.table(combined, glue("/data/xsun/tmp/cnv_{class}_g34.sum"), col.names = T, row.names = T, quote = F, sep = "\t")








