library(GenomicRanges)
library(rtracklayer)
library(maftools)
library(dplyr)
library(glue)

mutation_colors <- c(
  "dmr_Amp" = "#E41A1C",    # Red, Amplification
  "dmr_Del" = "#377EB8",    # Blue, Deletion
  "dmr_SNV/INDEL" = "#FF7F00", # Orange, SNV/INDEL
  "dmr_Amp; gene_Amp" = "#F9A6A1",   # Red, Gene Amplification
  "dmr_Del; gene_Del" = "#A5C8E4",   # Blue, Gene Deletion
  "gene_Amp" = "#F9E4E2",
  "gene_Del" = "#E4F2FC",
  "Missense_Mutation" = "#4DAF4A", # Green, Missense Mutation
  "Nonsense_Mutation" = "#F781BF", # Pink, Nonsense Mutation
  "Nonstop_Mutation" = "#A65628", # Brown, Nonstop Mutation
  "Splice_Site" = "#984EA3",  # Purple, Splice Site Mutation, 
  "Multi_Hit" = "grey25", 
  "Complex_Event" = "lightgrey",
  "Frame_Shift_Del" = "#E7298A", 
  "Frame_Shift_Ins" = "#66C2A5",
  "In_Frame_Del" = "#FFFF33",
  "In_Frame_Ins" = "#A6D854"
)

gistic_intersect <- function(gr, gistic_gr) {
  olp <- findOverlapPairs(gr, gistic_gr)
  res <- IRanges::pintersect(olp@first, olp@second)
  res$hit <- NULL
  return(res)
}

#' only return gr1 with over x percentage overlap either in gr1 or gr2
#' @param gr1 (GRanges)
#' @param gr2 (GRanges)
#' @return (GRanges)
perc_overlap <- function(gr1, gr2, thresh = 0.5) {
  olp <- findOverlapPairs(gr1, gr2)
  percentOverlap_gr1 <- width(IRanges::pintersect(olp@first, olp@second))/width(olp@first)
  percentOverlap_gr2 <- width(IRanges::pintersect(olp@first, olp@second))/width(olp@second)
  
  res <- olp@first[which(percentOverlap_gr1 > thresh | percentOverlap_gr2 > thresh)]
  
  return(res)
}

getGenes <- function(
    gene_types = NULL, 
    geneFile = "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"
    ) {
  genes <- readGFF(geneFile)
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
  return(geneGR)
}


#' get nearest gene to all ranges in gr
getNearestGene <- function(
    gr, 
    gene_types = NULL, 
    geneFile = "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"
    ) {
  genes <- readGFF(geneFile)
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

#' get G3/4 MB genes
#' @param db (character) The directory to MB_gene list
#' @return (character) A vector of Grp3/4 MB genes
get_g34genes <- function(db = glue("/.mounts/labs/pailab/src/gene_list/MB_gene",
                                   "/MBgene_database_20240709171001.csv")
) {
  df <- read.csv(db, stringsAsFactors = F, header = T, row.names = 1)
  sub <- df[,grepl("Hendrikse2022_G3_G4_MB_genes|Northcott2017_G34_genes", 
                   colnames(df)
  )
  ]
  
  g34_genes <- unique(rownames(sub[rowSums(sub) != 0,]))
  
  return(g34_genes)
}


#' get neurodevelopment genes
#' @param path (charaters) The path to the file containing the gene list
#' @return (character) A vector of neurodevelopment genes
get_neurodevGenes <- function() {
  additional_genes = c("OLIG3") # mentioned by Kim Aldinger
  
  paths <- c(avc = glue("/.mounts/labs/pailab/src/gene_list/brain-development",
                        "/aldinger2021_vladoiu2019_carter2018.csv"),
             leto = glue("/.mounts/labs/pailab/src/gene_list/brain-development",
                         "/LetoEtAl.txt"),
             bhaduri = glue("/.mounts/labs/pailab/src/gene_list/brain-development",
                            "/BhaduriKriegstein2021_nature_neocortex.txt"),
             sepp = glue("/.mounts/labs/pailab/src/gene_list/brain-development",
                         "/Sepp2023_CB.txt"),
             ian_leo = glue("/.mounts/labs/pailab/src/gene_list/brain-development",
                            "/cell_gene_mapping.csv")
  )
  
  # aldinger2021_vladoiu2019_carter2018.csv #
  avc <- read.csv(paths["avc"], stringsAsFactors = F, header = T)
  # select only stem/neuro progenitors
  avc_genes <- unique(avc[grepl("RL|UBCs|GCPs|stem", avc$region), "human"])
  
  # LetoEtAl.txt #
  leto <- read.table(paths["leto"], stringsAsFactors = F, header = F, sep = "\t")
  leto_genes <- unique(leto$V2)
  
  # BhaduriKriegstein2021_nature_neocortex.txt #
  bhaduri <- read.table(paths["bhaduri"], 
                        stringsAsFactors = F, header = F, sep = "\t"
  )
  bhaduri_genes <- bhaduri[! grepl("OPC|astrocytes", bhaduri$V2), "V1"]
  
  # Sepp2023 #
  sepp <- read.table(paths["sepp"], stringsAsFactors = F, header = F, sep = "\t")
  sepp_genes <- unlist(stringr::str_split(sepp$V2, ", "))
  sepp_genes <- stringr::str_trim(sepp_genes)
  sepp_genes <- unique(sepp_genes)
  
  # Ian & Leo #
  ianLeo <- read.csv(paths["ian_leo"], stringsAsFactors = F, header = T)
  ianLeo <- ianLeo[!ianLeo$reference %in% c("DEA", ""),]
  ianLeo <- ianLeo[!grepl("(oligodendrocytes|Microglia|Endothelial|Astrocytes|Purkinje)", 
                          ianLeo$lineage),]
  ianLeo_genes <- unique(ianLeo$gene)
  
  ### Combine ###
  res <- unique(c(avc_genes, leto_genes, bhaduri_genes, sepp_genes, ianLeo_genes,
                  additional_genes)
  )
  message(sprintf("Obtained %i neuro development genes from %s", 
                  length(res),
                  paste(c("", paths), collapse = "\n- ")
  )
  )
  
  return(res)
}


### gene ###
genes <- getGenes(gene_types = c("protein_coding"))

### g34/neurodev genes ###
g34_genes <- get_g34genes(glue("/.mounts/labs/pailab/src/gene_list/MB_gene",
                               "/MBgene_database_20240709171001.csv")
)
neurodev_genes <- get_neurodevGenes()

### muts ###
g4_muts <- import.bed("/.mounts/labs/pailab/private/xsun/output/ncMutMB/20240314/mutation/PCAWG/merged/Group4_PCAWG_snv-indel_hg38.bed")
g3_muts <- import.bed("/.mounts/labs/pailab/private/xsun/output/ncMutMB/20240314/mutation/PCAWG/merged/Group3_PCAWG_snv-indel_hg38.bed")

### gistic peaks ###
gistic_g4 <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/GenePattern/pcawg_group4_mb_gistic.regions_track.conf_95.hg38.bed")
gistic_g3 <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/GenePattern/pcawg_group3_mb_gistic.regions_track.conf_95.hg38.bed")

gistic_g4_amp <- gistic_g4[startsWith(gistic_g4$name, "Any-AP")]
gistic_g4_del <- gistic_g4[startsWith(gistic_g4$name, "Any-DP")]

gistic_g3_amp <- gistic_g3[startsWith(gistic_g3$name, "Any-AP")]
gistic_g3_del <- gistic_g3[startsWith(gistic_g3$name, "Any-DP")]


### cna ###
amp_g4 <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/pcawg_bed/Group4_pcawg.hg38.seg.amp.bed")
del_g4 <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/pcawg_bed/Group4_pcawg.hg38.seg.del.bed")

amp_g3 <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/pcawg_bed/Group3_pcawg.hg38.seg.amp.bed")
del_g3 <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/pcawg_bed/Group3_pcawg.hg38.seg.del.bed")

### dmr ###
# dmrs.bed was generated by reformatting the 
# /.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMRs.csv
# to bed using awk
dmr <- import.bed("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/dmrs.bed")
dmr <- getNearestGene(dmr, gene_types = c("protein_coding"))


#################################
summarize_onco <- function(dmr, amp, del, mut = NULL,
                           amp_gistic_gr_list = NULL,
                           del_gistic_gr_list = NULL,
                           chip = NULL, chip_ol_perc_thresh = 0.8,
                           genes = NULL
) {
  ### process dmr ###
  if (! is.null(chip)) {
    dmr <- perc_overlap(gr1 = dmr, gr2 = chip, thresh = chip_ol_perc_thresh)
  }
  
  ### process amp/del ###
  if (! is.null(amp_gistic_gr_list)) {
    for (gr in amp_gistic_gr_list) {
      amp_gistic <- gistic_intersect(amp, gr)
    }
  } else {
    amp_gistic <- amp
  }
  
  if (! is.null(del_gistic_gr_list)) {
    for (gr in del_gistic_gr_list) {
      del_gistic <- gistic_intersect(del, gr)
    }
  } else {
    del_gistic <- del
  }
  
  
  ### amp ###
  ## GSITIC dmr amp ##
  amp_olp <- findOverlapPairs(dmr, amp_gistic, type = "within")
  amp_dmr <- amp_olp@first
  amp_dmr$id <- amp_olp@second$name
  amp_dmr$CN <- "dmr_Amp"
  amp_dmr_df <- data.frame(Gene = amp_dmr$nearestGene, 
                           Sample_name = amp_dmr$id,
                           CN = amp_dmr$CN
  )
  
  ## genes related to GISTIC dmr amp ##
  amp_dmrGenes <- genes[genes$name %in% unique(amp_dmr_df$Gene)]
  amp_dmrGene_ol <- findOverlapPairs(amp_dmrGenes, amp, type = "within")
  amp_dmrGene <- amp_dmrGene_ol@first
  amp_dmrGene$id <- amp_dmrGene_ol@second$name
  amp_dmrGene$CN <- "gene_Amp"
  amp_dmrGene_df <- data.frame(Gene = amp_dmrGene$name, 
                               Sample_name = amp_dmrGene$id,
                               CN = amp_dmrGene$CN
  )
  
  ## GISTIC gene amp ##
  amp_gene_ol <- findOverlapPairs(genes, amp_gistic, type = "within")
  amp_gene <- amp_gene_ol@first
  amp_gene$id <- amp_gene_ol@second$name
  amp_gene$CN <- "gene_Amp"
  amp_gene_df <- data.frame(Gene = amp_gene$name, 
                            Sample_name = amp_gene$id,
                            CN = amp_gene$CN
  )
  
  amp_sum_df <- rbind(amp_dmr_df, amp_dmrGene_df, amp_gene_df)
  
  
  ### GISTIC dmr del ###
  del_olp <- findOverlapPairs(dmr, del_gistic, type = "within")
  del_dmr <- del_olp@first
  del_dmr$id <- del_olp@second$name
  del_dmr$CN = "dmr_Del"
  del_dmr_df <- data.frame(Gene = del_dmr$nearestGene, 
                           Sample_name = del_dmr$id,
                           CN = del_dmr$CN
  )
  
  ## genes related to GISTIC dmr del ##
  del_dmrGenes <- genes[genes$name %in% unique(del_dmr_df$Gene)]
  del_dmrGene_ol <- findOverlapPairs(del_dmrGenes, del)
  del_dmrGene <- del_dmrGene_ol@first
  del_dmrGene$id <- del_dmrGene_ol@second$name
  del_dmrGene$CN <- "gene_Del"
  del_dmrGene_df <- data.frame(Gene = del_dmrGene$name, 
                               Sample_name = del_dmrGene$id,
                               CN = del_dmrGene$CN
  )
  
  ## GISTIC gene del ##
  del_gene_ol <- findOverlapPairs(genes, del_gistic)
  del_gene <- del_gene_ol@first
  del_gene$id <- del_gene_ol@second$name
  del_gene$CN <- "gene_Del"
  del_gene_df <- data.frame(Gene = del_gene$name, 
                            Sample_name = del_gene$id,
                            CN = del_gene$CN
  )
  
  del_sum_df <- rbind(del_dmr_df, del_dmrGene_df, del_gene_df)
  
  
  ### muts ###
  if (! is.null(mut)) {
    mut_olp <- findOverlapPairs(dmr, mut)
    mut_sum <- mut_olp@first
    mut_sum$id <- mut_olp@second$name
    mut_sum$CN = "SNV/INDEL"
    mut_sum_df <- data.frame(Gene = mut_sum$nearestGene, 
                             Sample_name = mut_sum$id,
                             CN = "dmr_SNV/INDEL"
    )
    
    cna_sum_df <- rbind(amp_sum_df, del_sum_df, mut_sum_df)
  } else {
    cna_sum_df <- rbind(amp_sum_df, del_sum_df)
  }
  
  
  ### combine ###
  cna_sum_df <- 
    cna_sum_df %>% 
    unique() %>%
    group_by(Gene, Sample_name) %>%
    summarize(CN = paste(unique(CN), collapse = "; "), .groups = "drop")
  
  return(cna_sum_df)
}
#################################



  
### Grp 4###
g4_sum_df <- summarize_onco(
  dmr = dmr, amp = amp_g4, del = del_g4, mut = g4_muts,
  amp_gistic_gr_list = list(gistic_g4_amp, northcott_gistic_g4_amp),
  del_gistic_gr_list = list(gistic_g4_del, northcott_gistic_g4_del),
  #chip = g4_chip, chip_ol_perc_thresh = 0.9,
  genes = genes
  )


g4_amp_genes <- g4_sum_df %>% 
  filter(CN =="dmr_Amp") %>% 
  group_by(Gene) %>% 
  summarize(c = length(unique(Sample_name))) %>% 
  arrange(desc(c))

g4_del_genes <- g4_sum_df %>% 
  filter(CN =="dmr_Del") %>% 
  group_by(Gene) %>% 
  summarize(c = length(unique(Sample_name))) %>% 
  arrange(desc(c))

g4_dmr_mut <- g4_sum_df %>% 
  filter(CN =="dmr_SNV/INDEL") %>% 
  group_by(Gene) %>% 
  summarize(c = length(unique(Sample_name))) %>% 
  arrange(desc(c)) %>%
  filter(c > 1)


g4_maf_df_list <- lapply(
  list.files("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/g4_pcawg_maf", pattern = "*.maf$", full.names = T), 
  function(x) {
    tmp <- data.table::fread(x)
    #tmp$Tumor_Sample_Barcode <- stringr::str_split(x, "[/|.]", simplify = T)[,5]
    return(tmp)
    })

g4_onco <- read.maf(maf = data.table::rbindlist(g4_maf_df_list),
                        cnTable = g4_sum_df, 
                        verbose = FALSE)


oncoplot(maf = g4_onco, 
         genes = c(g34_genes[g34_genes %in% g4_onco@gene.summary$Hugo_Symbol]), 
         colors = mutation_colors
         )
oncoplot(maf = g4_onco, 
         genes = neurodev_genes[neurodev_genes %in% g4_onco@gene.summary$Hugo_Symbol],
         colors = mutation_colors
         )
oncoplot(maf = g4_onco, altered = T, 50, colors = mutation_colors)
oncoplot(maf = g4_onco, altered = T, 
         genes = c("CBFA2T3", "CBFA2T2", "PRDM6", "KDM6A", "RUNX1T1", "KDM2B"), 
         colors = mutation_colors)
oncoplot(maf = g4_onco, pathways = "sigpw", 
         gene_mar = 8, fontSize = 0.6, topPathways = 5, altered = T,
         colors = mutation_colors)
oncoplot(maf = g4_onco, altered = T, 
         genes = c("CBFA2T3", "CBFA2T2", "PRDM6", "KDM6A", "RUNX1T1", "KDM2B", "RPTOR"),
         colors = mutation_colors)
oncoplot(maf = g4_onco, 
         genes = c(g4_amp_genes$Gene, g4_del_genes$Gene, g4_dmr_mut$Gene), 
         colors = mutation_colors)

oncoplot(maf = g4_onco, pathways = "smgbp", 
         gene_mar = 8, fontSize = 0.6, topPathways = 5, altered = T,
         colors = mutation_colors)

mTORC1_pathway_genes <- c("AKT1","AKT1S1","ATG13","BNIP3","BRAF","CCNE1","CDK2","CLIP1","CYCS","DDIT4","DEPTOR","EEF2","EIF4A1","EIF4B","EIF4E","EIF4EBP1","FBXW11","HRAS","IKBKB","IRS1","MAP2K1","MAP2K2","MAPK1","MAPK3","MAPKAP1","MLST8","MTOR","NRAS","PDCD4","PDPK1","PLD1","PLD2","PML","POLDIP3","PPARGC1A","PRKCA","PRR5","PXN","RAC1","RAF1","RB1CC1","RHEB","RHOA","RICTOR","RPS6KA1","RPS6KB1","RPTOR","RRAGA","RRAGB","RRAGC","RRAGD","RRN3","SFN","SGK1","SREBF1","SSPO","TSC1","TSC2","ULK1","ULK2","YWHAB","YWHAE","YWHAG","YWHAH","YWHAQ","YWHAZ","YY1")
oncoplot(maf = g4_onco, 
         genes = c(mTORC1_pathway_genes[mTORC1_pathway_genes %in% g4_onco@gene.summary$Hugo_Symbol]), 
         colors = mutation_colors)


### Grp 3 ###
g3_sum_df <- summarize_onco(
  dmr = dmr, amp = amp_g3, del = del_g3, mut = g3_muts,
  amp_gistic_gr_list = list(gistic_g3_amp, northcott_gistic_g3_amp),
  del_gistic_gr_list = list(gistic_g3_del, northcott_gistic_g3_del),
  #chip = g3_chip, chip_ol_perc_thresh = 0.9,
  genes = genes
)


g3_amp_genes <- g3_sum_df %>% 
  filter(CN =="dmr_Amp") %>% 
  group_by(Gene) %>% 
  summarize(c = length(unique(Sample_name))) %>% 
  arrange(desc(c))

g3_del_genes <- g3_sum_df %>% 
  filter(CN =="dmr_Del") %>% 
  group_by(Gene) %>% 
  summarize(c = length(unique(Sample_name))) %>% 
  arrange(desc(c))

g3_dmr_mut <- g3_sum_df %>% 
  filter(CN =="dmr_SNV/INDEL") %>% 
  group_by(Gene) %>% 
  summarize(c = length(unique(Sample_name))) %>% 
  arrange(desc(c)) %>%
  filter(c > 1)


g3_maf_df_list <- lapply(
  list.files("/.mounts/labs/pailab/private/xsun/output/pcawgMAF/output/g3_pcawg_maf", pattern = "*.maf$", full.names = T), 
  function(x) {
    tmp <- data.table::fread(x)
    #tmp$Tumor_Sample_Barcode <- stringr::str_split(x, "[/|.]", simplify = T)[,5]
    return(tmp)
  })

g3_onco <- read.maf(maf = data.table::rbindlist(g3_maf_df_list),
                    cnTable = g3_sum_df, 
                    verbose = FALSE)
oncoplot(maf = g3_onco, 
         genes = c(mTORC1_pathway_genes[mTORC1_pathway_genes %in% g4_onco@gene.summary$Hugo_Symbol]), 
         colors = mutation_colors)

oncoplot(maf = g3_onco, 
         genes = c(g34_genes[g34_genes %in% g3_onco@gene.summary$Hugo_Symbol]), 
         colors = mutation_colors
)

oncoplot(maf = g3_onco, 
         genes = neurodev_genes[neurodev_genes %in% g3_onco@gene.summary$Hugo_Symbol],
         colors = mutation_colors
)

oncoplot(maf = g3_onco,
         genes = c("MYC","GFI1B","SMARCA4","KBTBD4","CTDNEP1","KMT2D","MYCN",
                   "GFI1","OTX2","BRCA2", 
                   
                   )
         )

oncoplot(maf = g3_onco, pathways = "sigpw", 
         gene_mar = 8, fontSize = 0.6, topPathways = 50, altered = T,
         colors = mutation_colors)

oncoplot(maf = g3_onco, pathways = "smgbp", 
         gene_mar = 8, fontSize = 0.6, topPathways = 5, altered = T,
         colors = mutation_colors)

oncoplot(maf = g3_onco, 
         genes = c(g3_amp_genes$Gene, g3_del_genes$Gene, g3_dmr_mut$Gene), 
         colors = mutation_colors)

oncoplot(maf = g3_onco, altered = T, 
         genes = c(g4_amp_genes$Gene, g4_del_genes$Gene, g4_dmr_mut$Gene),
         colors = mutation_colors)


### combined ###
coOncoplot(m1 = g4_onco, m2 = g3_onco,
           m1Name = "Grp4", m2Name = "Grp3",
           genes = g34_genes, 
           removeNonMutated = TRUE, 
           colors = mutation_colors, 
           sortByM1 = T
           )

coOncoplot(m1 = g4_onco, m2 = g3_onco,
           m1Name = "Grp4", m2Name = "Grp3",
           genes = neurodev_genes, 
           removeNonMutated = TRUE, 
           sortByM1 = T,
           colors = mutation_colors
)

coOncoplot(m1 = g4_onco, m2 = g3_onco,
           m1Name = "Grp4", m2Name = "Grp3",
           genes = c(g4_amp_genes$Gene, g4_del_genes$Gene, g4_dmr_mut$Gene,
                     g3_amp_genes$Gene, g3_del_genes$Gene, g3_dmr_mut$Gene), 
           removeNonMutated = TRUE, 
           colors = mutation_colors
)

coOncoplot(m1 = g4_onco, m2 = g3_onco,
           m1Name = "Grp4", m2Name = "Grp3",
           genes = c("CBFA2T3", "CBFA2T2", "PRDM6", "KDM6A", "RUNX1T1", "KDM2B", 
                     "KBTBD4", "GSE1", "IRF2BP2", "GFI1", "GFI1B"), 
           removeNonMutated = F, 
           colors = mutation_colors
)

coOncoplot(m1 = g4_onco, m2 = g3_onco,
           m1Name = "Grp4", m2Name = "Grp3",
           genes = mTORC1_pathway_genes, 
           removeNonMutated = TRUE, 
           colors = mutation_colors, sortByM1 = T
)

coOncoplot(m1 = g4_onco, m2 = g3_onco,
           m1Name = "Grp4", m2Name = "Grp3",
           genes = c("BARHL1", "DDX31", "MYC", "GFI1", "GFI1B"), 
           removeNonMutated = TRUE, 
           colors = mutation_colors, sortByM1 = T
)






## find top dmrs
#neurodev_genes <- read.table("input/neurodev.genes")[[1]]
#d <- c(amp_sum, del_sum)
#neurodev_genes <- read.table("input/neurodev.genes")[[1]]
#neurodev_genes <- read.table("input/neurodev.genes", header = F)[[1]]
#g34_genes <- read.table("input/g34mb.genes", header = F)[[1]]
#d$is_neurodev <- d$nearestGene %in% neurodev_genes
#d$is_mb <- d$nearestGene %in% g34_genes
#
#sub <- d[d$is_neurodev | d$is_mb]
#sub <- unique(sub)
#
#target_ranking <- read.table("./input/oligoDesign/targetRanking_240731.tsv", 
#                             stringsAsFactors = F, header = T, sep = "\t"
#                             )
#
#rank <- GenomicRanges::makeGRangesFromDataFrame(target_ranking, keep.extra.columns = T)
#tmp <- findOverlapPairs(sub, rank)
#res <- tmp@first
#res$rank <- tmp@second$rank
#write.table(as.data.frame(res), "~/Desktop/ranked_cna_dmrs.tsv", sep = "\t")
#
#
#