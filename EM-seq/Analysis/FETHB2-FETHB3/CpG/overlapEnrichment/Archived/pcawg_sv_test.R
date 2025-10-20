rm(list=ls())
library(BSgenome.Hsapiens.UCSC.hg38.masked) # needed for genNullSeqs
library(tidyr)
library(glue)
library(IRanges)


setwd(this.path::this.dir()) # set current scripts' dir as working dir
source("../../utils.R")
source("../../../../../MPRA/oligo_design/utils.R")
source("utils.R")

cpg_dmr_date <- get_configs("CPG_DMR_DATE")
projectRoot <- glue("/.mounts/labs/pailab/private/projects",
                    "/FetalHindbrain/EMseq_FETHB3/output/downstream",
                    "/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/",
                    "withoutBatchCorrection/{cpg_dmr_date}")


dt <- format(Sys.Date(),"%y%m%d")

dmr <- getDMRs(projectRoot)
dmr <- regioneR::filterChromosomes(dmr, organism="hg", chr.type="canonical")
cat(sprintf("Filter standard chrom: DMRs left=%i\n", length(dmr)))


tryCatch({
    # feCB enh
    fetal_CB_enh <- get_fetalCB_enh(sprintf("%s/raw",
                                            get_configs("ALDINGER_CB_CRE_DIR"))
    )
    
    # load tumour atac
    g3_atac <- import.bed(
      "/data/xsun/20240314/cre/Smith2022-Group3/allEnhancerLikeElements.bed"
    )
    g4_atac <- import.bed(
      "/data/xsun/20240314/cre/Smith2022-Group4/allEnhancerLikeElements.bed"
    )
    
    # load hg38 sv
    tmp <- read.table("/data/xsun/ol_test/wgs/pcawg_svs/hg38/Group34_sv_hg38.bed", 
                      sep = "\t", header = F)
    colnames(tmp) <- c("seqnames", "start", "end", "sample_id", "sv_class", "n_caller")
    tmp$start <- tmp$start + 1 # bed to 1-based GRanges
    pcawg_sv_hg38 <- GenomicRanges::makeGRangesFromDataFrame(tmp, keep.extra.columns = T)
    
    # add group info
    id <- read.csv("/data/xsun/ol_test/wgs/pcawg_svs/id_subgroup.csv")
    info <- id$Subgroup
    names(info) <- id$id
    pcawg_sv_hg38$subgroup <- info[pcawg_sv_hg38$sample_id]
    
    # split by svclass
    del <- pcawg_sv_hg38[pcawg_sv_hg38$sv_class == "DEL"]
    dup <- pcawg_sv_hg38[pcawg_sv_hg38$sv_class == "DUP"]
    inv <- pcawg_sv_hg38[pcawg_sv_hg38$sv_class %in% c("h2hINV", "t2tINV")]
    deldup <- pcawg_sv_hg38[pcawg_sv_hg38$sv_class %in% c("DEL", "DUP")]
    
    sapply(list(del, dup, inv, deldup), length)
    
    ### concurrent sequences
    g34_genes <- get_g34genes()
    tss <- get_gencode_anno()
    g34_p <- tss[tss$name %in% g34_genes]
    
    #gr <- pcawg_sv_hg38
    gr <- deldup
    #gr <- del
    print(length(gr))
    #gr <- gr[gr$n_caller >= 3]
    gr <- gr[width(gr) >= 100 & width(gr) <= 5e8]
    print(length(gr))
    
    
    # calc coverage of gr
    gr_split <- split(gr, gr$sample_id)
    
    gr_split <- lapply(gr_split, function(x) {
    condition <- length(findOverlaps(x[x$sv_class %in% c("DUP", "DEL")], g34_p)) > 0
    #condition <- F
    if (condition) {
      return(GRanges())
    } else {
      return(x)
    }
    })
    
    quantile(unlist(lapply(gr_split, length)))
    print(length(unlist(GRangesList(gr_split))))
    
    gr_split_reduced <- GRangesList(lapply(gr_split, reduce))
    quantile(unlist(lapply(gr_split_reduced, length)))
    
    cov <- IRanges::coverage(gr_split_reduced) # coverage by sample
    covered_ranges <- IRanges::slice(cov, lower=3, rangesOnly=T)
    res <- GRanges(covered_ranges)
    res
    
    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(res, dmr)), start(subsetByOverlaps(res, dmr)), end(subsetByOverlaps(res, dmr)))
    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(res, hars)), start(subsetByOverlaps(res, hars)), end(subsetByOverlaps(res, hars)))
    
    dmr_resized <- resize(dmr, width = width(dmr)+2e4*2, fix = "center")
    hars_resized <- resize(hars, width = width(hars)+2e4*2, fix = "center")
    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(res, dmr_resized)), start(subsetByOverlaps(res, dmr_resized)), end(subsetByOverlaps(res, dmr_resized)))
    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(dmr_resized, res)), start(subsetByOverlaps(dmr_resized, res)), end(subsetByOverlaps(dmr_resized, res)))
    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(res, hars_resized)), start(subsetByOverlaps(res, hars_resized)), end(subsetByOverlaps(res, hars_resized)))
    
    
    # summary event type of recurrent regions ol dmrs
    res_ol_dmr <- subsetByOverlaps(res, IRanges::intersect(dmr, c(g3_atac, g4_atac)))
    res_ol_dmr
    
    tmp <- res_ol_dmr[1]
    table(subsetByOverlaps(gr, tmp)$sv_class, subsetByOverlaps(gr, tmp)$sample_id)
    table(subsetByOverlaps(gr, tmp)$sv_class, subsetByOverlaps(gr, tmp)$subgroup)

    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(dmr, tmp)), start(subsetByOverlaps(dmr, tmp)), end(subsetByOverlaps(dmr, tmp)))
    
    
    # summary event type of recurrent regions ol hars
    res_ol_hars <- subsetByOverlaps(res, IRanges::intersect(hars, c(g3_atac, g4_atac)))
    res_ol_hars
    
    tmp <- res_ol_hars[1]
    table(subsetByOverlaps(gr, tmp)$sv_class, subsetByOverlaps(gr, tmp)$sample_id)
    table(subsetByOverlaps(gr, tmp)$sv_class, subsetByOverlaps(gr, tmp)$subgroup)
    
    sprintf("%s:%d-%d", seqnames(subsetByOverlaps(hars, tmp)), start(subsetByOverlaps(hars, tmp)), end(subsetByOverlaps(hars, tmp)))
    
    subsetByOverlaps(subsetByOverlaps(hars, tmp), g3_atac)
    subsetByOverlaps(subsetByOverlaps(hars, tmp), g4_atac)
    subsetByOverlaps(subsetByOverlaps(hars, tmp), fetal_CB_enh)
    
    
    ids <- subsetByOverlaps(gr[gr$sv_class == "DUP" & gr$subgroup == "Group4"], tmp)$sample_id
    id[id$id %in% ids,]
    gr_split[ids]
    
})
