rm(list=ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(glue)

default_wd <- this.path::this.dir()
setwd(default_wd) # set current scripts' dir as working dir
source("../../EM-seq/Analysis/FETHB2-FETHB3/utils.R")
source("../../EM-seq/Analysis/FETHB2-FETHB3/CpG/overlapEnrichment/utils.R")
source("./utils.R")

dt <- format(Sys.Date(),"%y%m%d")
dmr_date <- get_configs("CPG_DMR_DATE")
outDir <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/{dmr_date}/oligoDesign")

logFile <- sprintf("%s/target_ranking_%s.log",
                   outDir, 
                   dt
)

configs <- yaml::read_yaml("./config.yaml")

targets <- read.table(sprintf("%s/trimmed_targets_%s.tsv", outDir, "240731"))


main <- function() {
  targets <- GenomicRanges::makeGRangesFromDataFrame(targets, 
                                                     keep.extra.columns = T
  )
  
  ### Data Prep for annotation ###
  ## Whalen N2/N3 ##
  N2_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                          "/ChIP-seq/Whalen_2023",
                          "/GSE110758_human-HS1-11-N2-pooled-K27ac.narrowPeak.gz")
  N3_h3k27ac_file <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                          "/ChIP-seq/Whalen_2023",
                          "/GSE110758_human-HS1-11-N3-pooled-K27ac.narrowPeak.gz")
  N2_h3k27ac <- rtracklayer::import(N2_h3k27ac_file)
  N2_h3k27ac <- liftOver_gr(N2_h3k27ac)
  
  N3_h3k27ac <- rtracklayer::import(N3_h3k27ac_file)
  N3_h3k27ac <- liftOver_gr(N3_h3k27ac)
  
  # N2
  overlapping_targets <- perc_overlap(targets, N2_h3k27ac)
  targets_ol_N2h3k27ac <- paste(seqnames(overlapping_targets), 
                             start(overlapping_targets), 
                             end(overlapping_targets), 
                             sep = "_"
  )
  message(
    sprintf(
      "%d targets (duplication kept; %d unique targets) retrieved from N2_h3k27ac", 
      length(targets_ol_N2h3k27ac), 
      length(unique(targets_ol_N2h3k27ac))
    )
  )
  
  # N3
  overlapping_targets <- perc_overlap(targets, N3_h3k27ac)
  targets_ol_N3h3k27ac <- paste(seqnames(overlapping_targets), 
                                start(overlapping_targets), 
                                end(overlapping_targets), 
                                sep = "_"
  )
  message(
    sprintf(
      "%d targets (duplication kept; %d unique targets) retrieved from N3_h3k27ac", 
      length(targets_ol_N3h3k27ac), 
      length(unique(targets_ol_N3h3k27ac))
    )
  )
  
  
  ## dmr ol Aldinger human fetal hindbrain histone modifications (ChIP-seq) ##
  fetal_CB_enh <- get_fetalCB_enh(sprintf("%s/raw",
                                          get_configs("ALDINGER_CB_CRE_DIR"))
                                  )
  
  ## fetal hindbrain enhancer overlap ##
  overlapping_targets <- perc_overlap(targets, fetal_CB_enh)
  targets_ol_fcbEnh <- paste(seqnames(overlapping_targets), 
                             start(overlapping_targets), 
                             end(overlapping_targets), 
                             sep = "_"
  )
  message(
    sprintf(
      "%d targets (duplication kept; %d unique targets) retrieved from Aldinger fetal hindbrain enhancer overlap", 
      length(targets_ol_fcbEnh), 
      length(unique(targets_ol_fcbEnh))
    )
  )
  
  
  
  ## SV ##
  sv <- getNorthcott2012_AmpsDels()
  targets_ol_G3ampSV <- get_ol_targets(targets, 
                                       sv$amps$`GISTIC_Amps-Group3`, 
                                       "G3 MB ampSV"
  )
  targets_ol_G4ampSV <- get_ol_targets(targets,
                                       sv$amps$GISTIC_Amps_Group4,
                                       "G4 MB ampSV"
  )
  targets_ol_G3delSV <- get_ol_targets(targets,
                                       sv$dels$GISTIC_Dels_Group3, 
                                       "G3 MB delSV"
  )
  targets_ol_G4delSV <- get_ol_targets(targets,
                                       sv$dels$GISTIC_Dels_Group4, 
                                       "G4 MB delSV"
  )
  
  ## HAR ##
  har <- getHARs()
  targets_ol_hars <- get_ol_targets(targets, har, "HARs")
  
  ## ENCODE enhancer ##
  encode_ccres <- get_encode_cres()
  encode_ELS <- encode_ccres[grepl("ELS", encode_ccres$Type)]
  targets_ol_encodeELS <- get_ol_targets(targets, 
                                         encode_ccres, 
                                         "ENCODE Enhancer-like elements"
  )
  
  ## known G34 MB genes ##
  g34_genes <- get_g34genes(glue("/.mounts/labs/pailab/src/gene_list/MB_gene",
                                 "/MBgene_database_20240709171001.csv")
                            )
  
  ## G34 SNV/INDEL (PEMECA-PCAWG) ##
  g34_mut <- import.bed(glue("/data/xsun/20240314/mutation/PEMECA-PCAWG/merged",
                             "/Group34_PEMECA-PCAWG_snv-indel_hg38.bed"
                             )
                        )
  overlapping_targets <- perc_overlap(targets, g34_mut)
  targets_ol_mut <- paste(seqnames(overlapping_targets), 
                          start(overlapping_targets), 
                          end(overlapping_targets), 
                          sep = "_"
  )
  message(
    sprintf(
      "%d targets (duplication kept; %d unique targets) retrieved from G34 mutation (PEMECA-PCAWG) overlap", 
      length(targets_ol_mut), 
      length(unique(targets_ol_mut))
    )
  )
  
  ## Smith et al MB G3 atac ##
  g3_atac <- import.bed("/data/xsun/20240314/cre/Smith2022-Group3/allEnhancerLikeElements.bed")
  
  overlapping_targets <- perc_overlap(targets, g3_atac)
  targets_ol_g3atac <- paste(seqnames(overlapping_targets), 
                             start(overlapping_targets), 
                             end(overlapping_targets), 
                             sep = "_"
  )
  message(
    sprintf(
      "%d targets (duplication kept; %d unique targets) retrieved from Smith et al G3 ATAC overlap", 
      length(targets_ol_g3atac), 
      length(unique(targets_ol_g3atac))
    )
  )
  
  ## Smith et al MB G4 atac ##
  g4_atac <- import.bed("/data/xsun/20240314/cre/Smith2022-Group4/allEnhancerLikeElements.bed")
  
  overlapping_targets <- perc_overlap(targets, g4_atac)
  targets_ol_g4atac <- paste(seqnames(overlapping_targets), 
                             start(overlapping_targets), 
                             end(overlapping_targets), 
                             sep = "_")
  message(
    sprintf(
      "%d targets (duplication kept; %d unique targets) retrieved from Smith et al G4 ATAC overlap", 
      length(targets_ol_g4atac), 
      length(unique(targets_ol_g4atac))
    )
  )
  
  
  ## DepMap genes ##
  # DepMap G3 MB cell lines: D341Med, D283MED, D425, D458
  # Downloaded from DepMap portal: https://depmap.org/portal/data_page/?tab=customDownloads
  depmap_g3_mb <- read.csv(
    glue("/data/xsun/db/DepMap",
         "/CRISPR_(DepMap_Public_24Q2+Score,_Chronos)_subsetted_NAsdropped.csv"
    ), 
    stringsAsFactors = F, 
    header = T
  )
  median_score <- apply(depmap_g3_mb[, 9:ncol(depmap_g3_mb)], 2, median)
  
  # neuroDev genes #
  neurodev_genes <- get_neurodevGenes()
  
  ## ABC prediction ##
  abc <- get_neuroABC()
  ol <- findOverlaps(targets, abc)
  overlapping_targets <- targets[queryHits(ol)]
  overlapping_abc <- abc[subjectHits(ol)]
  
  targets_ol_abc_df <- as.data.frame(overlapping_targets)
  targets_ol_abc_df$target <- paste(targets_ol_abc_df$seqnames, 
                                    targets_ol_abc_df$start,
                                    targets_ol_abc_df$end,
                                    sep = "_"
  )
  targets_ol_abc_df$abcTarget <- overlapping_abc$name
  targets_ol_abc_df$depmapScore <- median_score[targets_ol_abc_df$abcTarget]
  targets_ol_abc_df$isG34gene <- targets_ol_abc_df$abcTarget %in% g34_genes
  targets_ol_abc_df$isNeurodev <- targets_ol_abc_df$abcTarget %in% neurodev_genes
  targets_abcTargets <- targets_ol_abc_df %>% 
    dplyr::select(target, abcTarget, depmapScore, isG34gene, isNeurodev) %>% 
    group_by(target) %>%
    summarise(abcTargetGenes = toString(sort(unique(abcTarget))), 
              lowestAbcDepmapScore = ifelse(all(is.na(depmapScore)), 
                                            1e2, 
                                            min(depmapScore, na.rm = TRUE)
              ),
              abcTargetGenesIncludeMbGene = any(isG34gene),
              abcTargetGenesIncludeNeurodevGene = any(isNeurodev)
    )
  
  ## Geller2024 CRISPRi ##
  #' get Geller 2024 CRISPRi H9-derived neural stem cell growth assay results of conserved regions
  get_geller2024_crispri <- function() {
    library(liftOver)
    hg19_to_hg38_chainFile <- glue("/.mounts/labs/pailab/src/ucsc-tools",
                                   "/chain_files/hg19ToHg38.over.chain")
    
    conserved <- readxl::read_excel(glue("/.mounts/labs/pailab/src",
                                         "/neurodev-genomics/CRISPRiGrowthAssay",
                                         "/Geller_2024/mmc4.xlsx"), 
                                    sheet = 3
    )
    
    # Summarize status across the three time points #
    conserved <- conserved %>% mutate(pheno = case_when(
      grepl("dynamic", paste(t4.phenotype, t8.phenotype, t12.phenotype)) ~ "dynamic",
      grepl("positive", paste(t4.phenotype, t8.phenotype, t12.phenotype)) ~ "positive",
      grepl("negative", paste(t4.phenotype, t8.phenotype, t12.phenotype)) ~ "negative",
      .default = "neutral")
    )
    conserved$betaAbs <- apply(conserved[, c("t4.beta", "t8.beta", "t12.beta")], 
                               1, 
                               function(x) {max(abs(x))}
    )
    conserved <- conserved[, c("seqnames", "start", "end", "name", "sgRNA", "pheno", "betaAbs")]
    
    conserved_hg19_gr <- GenomicRanges::makeGRangesFromDataFrame(conserved, 
                                                                 keep.extra.columns = T
    )
    
    conserved_hg38_grList <- liftOver(conserved_hg19_gr, 
                                      import.chain(hg19_to_hg38_chainFile)
    )
    
    ln <- unlist(lapply(conserved_hg38_grList, length))
    if (any(ln!=1)) {
      # remove HARs with no perfect liftOver mapping
      conserved_hg38_grList <- conserved_hg38_grList[-which(ln!=1)] 
    }
    
    conserved_hg38_gr <- unlist(conserved_hg38_grList)
    
    cat(sprintf("%i out of %i converted to hg38\n", 
                length(conserved_hg38_gr), 
                length(conserved_hg19_gr)
    )
    )
    
    return(conserved_hg38_gr)
  }
  
  geller2024 <- get_geller2024_crispri()
  ol <- findOverlaps(targets, geller2024)
  overlapping_targets <- targets[queryHits(ol)]
  overlapping_geller <- geller2024[subjectHits(ol)]
  
  targets_ol_geller_df <- as.data.frame(overlapping_targets)
  targets_ol_geller_df$target <- paste(targets_ol_geller_df$seqnames, 
                                       targets_ol_geller_df$start,
                                       targets_ol_geller_df$end,
                                       sep
                                       = "_"
  )
  targets_ol_geller_df$geller_pheno <- overlapping_geller$pheno
  targets_ol_geller_df$geller_beta <- overlapping_geller$betaAbs
  
  targets_geller <- targets_ol_geller_df %>%
    dplyr::select(target, geller_pheno, geller_beta) %>% 
    group_by(target) %>%
    summarise(geller_crispri_pheno = toString(sort(unique(geller_pheno))), 
              geller_crispri_beta = max(geller_beta),
              geller_crispri_targets = n()
    )
  
  
  ### Annotate ###
  message("Annotating target regions...")
  # get nearestGene at the same time
  message("--- nearest gene")
  targets_df <- as.data.frame(
    getNearestGene(targets, gene_types = c("protein_coding"))
    ) 
  
  message("--- nearest gene DepMap score")
  targets_df$nearestGeneDepmapScore <- median_score[targets_df$nearestGene]
  # set NA with default 1e2
  targets_df$nearestGeneDepmapScore[is.na(targets_df$nearestGeneDepmapScore)] <- 1e2 
  targets_df$nearestGeneIsMbGene <- targets_df$nearestGene %in% g34_genes
  targets_df$nearestGeneIsNeurodevGene <- targets_df$nearestGene %in% neurodev_genes
  
  # get ABC targets info
  message("--- ABC neuronal cells target genes")
  targets_df <- merge(targets_df, targets_abcTargets, by = "target", all.x = T) 
  message("--- ABC neuronal cells target genes DepMap scores")
  targets_df$lowestAbcDepmapScore[is.na(targets_df$lowestAbcDepmapScore)] <- 1e2
  targets_df <- targets_df %>% 
    mutate(lowestDepmapScore = pmin(nearestGeneDepmapScore, 
                                    lowestAbcDepmapScore, 
                                    na.rm = T
    )
    )
  
  # get Geller CRISPRi info
  message("--- Geller H9-derived NSC CRISPRi beta")
  targets_df <- merge(targets_df, targets_geller, by = "target", all.x = T) 
  
  # not using count to avoid hypermutation
  message("--- PEMECA & PCAWG SNV/INDEL")
  targets_df$ol_SnvIndel <- table(targets_ol_mut)[targets_df$target]
  
  message("--- Northcott SV (amp/del)")
  targets_df$ol_G3ampSV <- table(targets_ol_G3ampSV)[targets_df$target]
  targets_df$ol_G4ampSV <- table(targets_ol_G4ampSV)[targets_df$target]
  targets_df$ol_G3delSV <- table(targets_ol_G3delSV)[targets_df$target]
  targets_df$ol_G4delSV <- table(targets_ol_G4delSV)[targets_df$target]
  targets_df$ol_G34SV <- table(
    c(targets_ol_G3ampSV, 
      targets_ol_G4ampSV, 
      targets_ol_G3delSV, 
      targets_ol_G4delSV)
  )[targets_df$target]

  message("--- HARs")
  targets_df$ol_har <- table(targets_ol_hars)[targets_df$target]
  
  message("--- ENCODE ELS")
  targets_df$ol_encodeELS <- targets_df$target %in% targets_ol_encodeELS
  
  # ol Aldinger fetal CB ChIP inferred enhancer
  message("--- Human fetal cerebellar enhancers inferred from Aldinger ChIP-seq")
  targets_df$ol_fcbEnh <- targets_df$target %in% targets_ol_fcbEnh
  
  # ol Whalen N2/N3 H3K27ac
  message("--- Whalen N2/N3 H3K27ac")
  targets_df$ol_N2h3k27ac <- targets_df$target %in% targets_ol_N2h3k27ac
  targets_df$ol_N3h3k27ac <- targets_df$target %in% targets_ol_N3h3k27ac
  
  # ol Smith G3 MB atac count
  message("--- Smith et al 2022 paper G3/4 MB ATAC peaks")
  targets_df$ol_g3atac <- targets_df$target %in% targets_ol_g3atac
  # ol Smith G4 MB atac count
  targets_df$ol_g4atac <- targets_df$target %in% targets_ol_g4atac
  
  count_cols <- c("ol_SnvIndel", "ol_G34SV",
                  "ol_G3ampSV", "ol_G4ampSV", "ol_G3delSV", "ol_G4delSV",
                  "ol_har"
                  
  )
  targets_df <- targets_df %>%
    mutate_at(vars(one_of(count_cols)), 
              ~ replace(., is.na(.), 0)) # fill NA counts by 0
  
  # promoter
  targets_notOL_promoters <- drop_promoters(
    GenomicRanges::makeGRangesFromDataFrame(targets_df, keep.extra.columns = T),
    method = "overlap",
    promoter_radius = 1000
  )$target
  
  targets_df$ol_promoter <- ! targets_df$target %in% targets_notOL_promoters
  
  
  ### Score and rank targets ###
  ## Score ##
  message("Ranking targets")
  score_dict <- c(
    ol_N2h3k27ac = 20,
    ol_N3h3k27ac = 20,
    ol_SnvIndel = 10, # the muts can be tested in MPRA
    ol_G3ampSV = 10, ol_G4ampSV = 10, ol_G3delSV = 10, ol_G4delSV = 10, 
    ol_har = 15, # we are interested in har given svz human evolution
    ol_encodeELS = 5, 
    ol_fcbEnh = 100, # for pilot we want them to be all fcbEnh
    ol_g3atac = 5, ol_g4atac = 5, # may drive tumour & consistent
    abcTargetGenesIncludeMbGene = 20, nearestGeneIsMbGene = 40, 
    abcTargetGenesIncludeNeurodevGene = 20, nearestGeneIsNeurodevGene = 40,
    geller_crispri_beta = 0
    )
  
  message(
    cat("Using score matrix as below:\n", 
        paste(names(score_dict), score_dict, "\n")
    )
  )
  
  score_df <- targets_df[, names(score_dict)]
  score_df <- apply(score_df, 2, as.numeric)
  score_df[is.na(score_df)] <- 0
  
  score_df_scaled <- apply(score_df, 2, function(x) {
    (x-min(x, na.rm = T))/((max(x, na.rm = T) - min(x, na.rm = T)))
  })
  
  score <- score_df_scaled %*% as.matrix(score_dict)
  
  ## Rank ##
  targets_df$score <- score[,1]
  targets_df <- targets_df %>%
    arrange(
      ol_promoter, # not ol promter first
      desc(score)
    ) %>%
    mutate(rank = row_number())
  
  outFile <- sprintf("%s/targetRanking_%s.tsv", outDir, dt)
  write.table(targets_df, outFile, col.names = T, row.names = T, quote = F, sep = "\t")
}



### main ###
if (! dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

logFileCon <- file(logFile, open = "wt")
sink(logFileCon, split = T, type = "output")
sink(logFileCon, type = "message")
tryCatch({main()}, 
         error = function(e) {message(e)},
         finally = {
           message("\n\n--------- R sessionInfo ---------\n\n")
           print(sessionInfo())
           sink(type = "output")
           sink(type = "message")
           close(logFileCon)
         }
) 