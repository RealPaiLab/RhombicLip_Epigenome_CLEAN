rm(list=ls())
library(UpSetR)


source("../../EM-seq/Analysis/FETHB2-FETHB3/utils_PaiLab.R")
source("../../EM-seq/Analysis/FETHB2-FETHB3/CpG/overlapEnrichment/utils_PaiLab.R")
source("./utils_PaiLab.R")

dt <- format(Sys.Date(),"%y%m%d")
dmr_date <- get_configs("CPG_DMR_DATE")
outDir <- glue("/home/rstudio/isilon/private/projects/FetalHindbrain",
               "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs",
               "/CTsnv_excluded/withoutBatchCorrection/240711/oligoDesign/crsAnnotate/250605")


configs <- yaml::read_yaml("./config.yaml")
mpra_target_len <- configs$MPRA_TARGET_LEN


### Load target crs ###
crsFile <- glue("/home/rstudio/isilon/private/projects/FetalHindbrain",
                "/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3",
                "/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711",
                "/oligoDesign/crsAnnotate/250605/targetProbes_hg38_250605.rds")

crs <- readRDS(crsFile)
class <- crs$class
source <- crs$source
elementMetadata(crs) <- NULL
crs$class <- class
crs$source <- source

browser()




### Whalen N2/N3 ###
N2_h3k27ac_file <- glue("/home/rstudio/isilon/src/neurodev-genomics",
                        "/ChIP-seq/Whalen_2023",
                        "/GSE110758_human-HS1-11-N2-pooled-K27ac.narrowPeak.gz")
N3_h3k27ac_file <- glue("/home/rstudio/isilon/src/neurodev-genomics",
                        "/ChIP-seq/Whalen_2023",
                        "/GSE110758_human-HS1-11-N3-pooled-K27ac.narrowPeak.gz")
N2_h3k27ac <- rtracklayer::import(N2_h3k27ac_file)
N2_h3k27ac <- liftOver_gr(N2_h3k27ac)

N3_h3k27ac <- rtracklayer::import(N3_h3k27ac_file)
N3_h3k27ac <- liftOver_gr(N3_h3k27ac)

## ol N2 ##
overlapping_crs <- perc_overlap(crs, N2_h3k27ac)
crs_ol_N2h3k27ac <- paste(seqnames(overlapping_crs), 
                              start(overlapping_crs), 
                              end(overlapping_crs), 
                              sep = "_"
)
message(
  sprintf(
    "%d crs (duplication kept; %d unique crs) retrieved from N2_h3k27ac", 
    length(crs_ol_N2h3k27ac), 
    length(unique(crs_ol_N2h3k27ac))
  )
)

## ol N3 ##
overlapping_crs <- perc_overlap(crs, N3_h3k27ac)
crs_ol_N3h3k27ac <- paste(seqnames(overlapping_crs), 
                              start(overlapping_crs), 
                              end(overlapping_crs), 
                              sep = "_"
)
message(
  sprintf(
    "%d crs (duplication kept; %d unique crs) retrieved from N3_h3k27ac", 
    length(crs_ol_N3h3k27ac), 
    length(unique(crs_ol_N3h3k27ac))
  )
)


### Nearest gene ###
crs <- getNearestGene(gr = crs, gene_types = "protein_coding")

## isNeurodev ##
isNeurodev <- crs$nearestGene %in% get_neurodevGenes()

## isG34MB ##
isG34MB <- crs$nearestGene %in% get_g34genes()


### ol feCB enh ###
feCBenh <- get_fetalCB_enh(sprintf("%s/raw",
                                        get_configs("ALDINGER_CB_CRE_DIR"))
                                )

overlapping_crs <- perc_overlap(crs, feCBenh)
crs_ol_feCBenh <- paste(seqnames(overlapping_crs), 
                          start(overlapping_crs), 
                          end(overlapping_crs), 
                          sep = "_"
                          )

### ol har ###
hars <- getHARs()
overlapping_crs <- perc_overlap(crs, hars)
crs_ol_hars <- paste(seqnames(overlapping_crs), 
                        start(overlapping_crs), 
                        end(overlapping_crs), 
                        sep = "_"
                        )

### ol K27ac summit ###
summit <- get_fetalCB_h3k27ac_summit(
  up = floor(mpra_target_len/2), 
  down = mpra_target_len - floor(mpra_target_len/2)
)
overlapping_crs <- perc_overlap(crs, summit)
crs_ol_summit <- paste(seqnames(overlapping_crs), 
                     start(overlapping_crs), 
                     end(overlapping_crs), 
                     sep = "_"
                     )

### SV ###
sv <- getNorthcott2012_AmpsDels()
g34_sv <- c(sv$amps$`GISTIC_Amps-Group3`, sv$amps$`GISTIC_Amps-Group4`,
            sv$dels$GISTIC_Dels_Group3, sv$dels$GISTIC_Dels_Group4)
overlapping_crs <- perc_overlap(crs, g34_sv)
crs_ol_sv <- paste(seqnames(overlapping_crs), 
                       start(overlapping_crs), 
                       end(overlapping_crs), 
                       sep = "_")

### SNV ###
g34_mut <- import.bed(glue("/home/rstudio/isilon/private/projects/FetalHindbrain/anno/Group34_PEMECA-PCAWG_snv-indel_hg38.bed"))
overlapping_crs <- perc_overlap(crs, g34_mut)
crs_ol_mut <- paste(seqnames(overlapping_crs), 
                        start(overlapping_crs), 
                        end(overlapping_crs), 
                        sep = "_"
                    )

### feCB H3K27ac UNION ###
fCB_h3k27ac <- getFetalCB_HistonePeaks()$H3K27ac_union
overlapping_crs <- perc_overlap(crs, fCB_h3k27ac)
crs_ol_feCBh3k27ac <- paste(seqnames(overlapping_crs), 
                    start(overlapping_crs), 
                    end(overlapping_crs), 
                    sep = "_"
)


### Sum & Plot ###
# class: N2/N3 K27ac; isNeurodev; isG34; feCB enh; HAR, K27ac summit
crs$id <- paste(seqnames(crs), start(crs), end(crs), sep = "_")

crs$ol_feCBenh <- as.numeric(crs$id %in% crs_ol_feCBenh)
crs$ol_feCBk27ac <- as.numeric(crs$id %in% crs_ol_feCBh3k27ac)
crs$ol_feCBk27ac_summit <- as.numeric(crs$id %in% crs_ol_summit)
crs$ol_N2_h3k27ac <- as.numeric(crs$id %in% crs_ol_N2h3k27ac)
crs$ol_N3_h3k27ac <- as.numeric(crs$id %in% crs_ol_N3h3k27ac)
crs$ol_hars <- as.numeric(crs$id %in% crs_ol_hars)
crs$isNeurodev <- as.numeric(isNeurodev)
crs$isG34MB <- as.numeric(isG34MB)
crs$ol_G34_sv <- as.numeric(crs$id %in% crs_ol_sv)
crs$ol_G34_snvIndel <- as.numeric(crs$id %in% crs_ol_mut)

crs$Neurodev_G34MB <- as.numeric(isNeurodev|isG34MB)
crs$N2_N3_k27ac <- as.numeric(crs$id %in% c(crs_ol_N2h3k27ac, crs_ol_N3h3k27ac))
crs$G34_mutated <- as.numeric(crs$id %in% c(crs_ol_sv, crs_ol_mut))


outFile <- sprintf("%s/oligo_summary_%s.tsv", outDir, dt)
write.table(as.data.frame(crs), outFile, 
            col.names = T, row.names = T, quote = F, sep = "\t")


intersted_cols <- c("ol_feCBenh", "ol_feCBk27ac", "ol_feCBk27ac_summit",
                    "N2_N3_k27ac",
                    "Neurodev_G34MB", 
                    "G34_mutated",
                    "ol_hars")

outFile <- sprintf("%s/oligo_summary_upSet_%s.png", outDir, dt)
png(outFile,width=10,height=6,units="in",res=600, pointsize=6)
upset(as.data.frame(crs), 
      sets = intersted_cols,
      order.by = "freq",
      text.scale = 1.6, point.size = 3,
      sets.bar.color = c("#f8101b", "olivedrab", "#f8101b","#fece2f",
                         "skyblue", "grey25", "#8e529e")
      )
dev.off()


outFile <- sprintf("%s/oligo_summary_bar_%s.png", outDir, dt)
png(outFile,width=7,height=4,units="in",res=600, pointsize=6)
barplot(sort(colSums(as.data.frame(crs)[, intersted_cols])/length(crs)*100, 
             decreasing = T), 
        ylab = "Percent CRS (%)", 
        ylim = c(0, 100)
        )
dev.off()


