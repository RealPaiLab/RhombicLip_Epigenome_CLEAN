library(glue)
library(GenomicRanges)

#### Utils for DMR analyses ####

#' get preset configs 
#' @param config_name (char) Name of the defined config.
get_configs <- function(config_name = NA) {
  ## CONFIGS ##
  # Cytosine reports #
  CPG_REPORT_DIR <- c(
    glue("/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3",
         "/output/alignment/methyldackel/report/CpG_snpFiltered")
  )
  
  # DMR #
  #CPG_DMR_DATE <- "240522" # outdated
  CPG_DMR_DATE <- "240711"
  CPG_DMR_FILE <- glue(
    "/.mounts/labs/pailab/private/projects/FetalHindbrain/EMseq_FETHB3",
    "/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded",
    "/withoutBatchCorrection/{CPG_DMR_DATE}/DMRs.csv"
  )
  
  # cCRE #
  ENCODE_CRE_FILE <- glue("/.mounts/labs/pailab/private/projects",
                          "/FetalHindbrain/anno/GRCh38-cCREs.bed")
  ALDINGER_CB_CRE_DIR <- glue("/.mounts/labs/pailab/private/xsun/output",
                              "/ncMutMB/20240314/cre/Aldinger-FetalCB")
  
  # Gene annotation #
  GENCODE_GENE_FILE <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                            "/anno/gencode.v44.basic.annotation.gtf")
  
  # Transposable elements #
  REPEATMASKER_TE_FILE <- glue("/.mounts/labs/pailab/private/xsun/Database",
                               "/RepeatMasker/RepeatMasker_open406_Dec2013_Dfam20_hg38",
                               "/hg38.fa.out.gz")
  
  # metadata #
  METADATA_FILE <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                        "/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv")
  
  # rl meme db #
  RL_MEME_DB <- glue("/data/xsun/db/meme/HOCOMOCOv12/filteredTFdb",
                     "/Hendrikse2022_RL_activeGenes_TF.meme")
  TF_DB <- "HOCOMOCOv12"
  
  # HOCOMOCOv12 
  HOCOMOCO_V12 <- glue("/.mounts/labs/pailab/src/transcription-factors",
                       "/HOCOMOCOv12/H12CORE_meme_format.meme")
  
  # Aldinger LCM RL-VZ vs SVZ RNA logFC
  RL_DEG <- glue("/.mounts/labs/pailab/private/projects/FetalHindbrain",
                "/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/231123",
                "/edgeR_RLVZvsSVZ_231123.txt")
  
  # Neuronal cell ABC targets
  NEURO_ABC_TARGETS <- glue("/data/xsun/20240314/creTargets/Nasser-Neuronal-ABC",
                            "/Nasser-Neuronal-ABC_creTarget_hg38.bed")
  
  # Whalen MPRA results
  WHALEN_MPRA_DIR <- glue("/.mounts/labs/pailab/src/neurodev-genomics/MPRA", 
                          "/Whalen_2023/GSE110758_RAW")
  
  # RL active genes
  RL_ACTIVE_GENES <- glue("/data/xsun/db/meme/HOCOMOCOv12/activeGenes",
                          "/Hendrikse2022_RL_activeGenes")
  
  # JASPAR
  JASPAR_DB <- "/.mounts/labs/pailab/src/ucsc-goldenpath/hg38/JASPAR2024.bb"
  
  # VISTA positive human hindbrain enhancer
  VISTA_HB <- glue("/.mounts/labs/pailab/src/neurodev-genomics",
                   "/TransgenicMouseAssay/VISTA_Enhancer_Browser",
                   "/hindbrain_positive_enhancers_hg19.vista")
  
  ## LIST ##
  configs = list(
    "CPG_REPORT_DIR"= CPG_REPORT_DIR, 
    "CPG_DMR_DATE" = CPG_DMR_DATE,
    "CPG_DMR_FILE" = CPG_DMR_FILE,
    "ENCODE_CRE_FILE" = ENCODE_CRE_FILE,
    "ALDINGER_CB_CRE_DIR" = ALDINGER_CB_CRE_DIR,
    "GENCODE_GENE_FILE" = GENCODE_GENE_FILE,
    "REPEATMASKER_TE_FILE" = REPEATMASKER_TE_FILE,
    "METADATA_FILE" = METADATA_FILE,
    "RL_MEME_DB" = RL_MEME_DB, 
    "HOCOMOCO_V12" = HOCOMOCO_V12,
    "TF_DB" = TF_DB,
    "RL_DEG" = RL_DEG,
    "NEURO_ABC_TARGETS" = NEURO_ABC_TARGETS, 
    "WHALEN_MPRA_DIR" = WHALEN_MPRA_DIR,
    "RL_ACTIVE_GENES" = RL_ACTIVE_GENES, 
    "JASPAR_DB" = JASPAR_DB,
    "VISTA_HB" = VISTA_HB
    )
  
  ## OUTPUT ##
  if (config_name %in% names(configs)) {
    return(configs[[config_name]])
  } else {
    config_names_str <- paste({names(configs)}, collapse = ", ")
    stop(glue(
      "Valid config_name is one of {config_names_str}."
    ))
  }
  
}


#' get EM-seq sample metadata
#' @return (data.frame) The metadata dataframe
get_sample_meta <- function() {
  meta <- read.table(get_configs("METADATA_FILE"), 
                     stringsAsFactors = F, 
                     header = T, 
                     sep = "\t"
                     )
  rownames(meta) <- meta$library_ID
  
  return(meta)
}


#' Extract EM-seq sample id from the given string
#' @param string (char) The target string.
#' @return (char) The sample id
extract_sampleId <- function(string) {
    return(stringr::str_extract(string, 
                                pattern = "FETHB[0-9]_[0-9]{4}_[0-9]{2}_LB0[12]"
                                )
           )
}


#' Get tissue type of the given sample id
#' @param sampleId (char) The sample id(s)
#' @return (numeric) Numeric label of VZ or SVZ tissue type
get_TissueType <- function(sampleId) {
  # VZ = 1; SVZ = 0
  tissues <- get_sample_meta()[sampleId, "ROI"]
  
  tissues[tissues == "VZ"] <- 1
  tissues[tissues == "SVZ"] <- 0
  tissues <- as.numeric(tissues)
  
  names(tissues) <- sampleId
  
  return(tissues)
}


#' Read in cytosine reports under provided directory/directories using bsseq::read.bismark
#' @param inDir (character) The directory/directories containing target cytosine reports
#' @param grepPattern (character) Regex to help identify the cytosine reports
#' @param excludeSamp (character) Sample names to be excluded (! not the file name).
#' @param drop_hg38_alt (logical) Remove alternative chromosomes in the hg38 genome build.
#' @return (list) Return a list containing filenames, sample names, and a bsseq object
readBS <- function(inDir, 
                   grepPattern="", 
                   excludeSamp=c(), 
                   drop_hg38_alt = TRUE){
	### Preprocess files ###
	fNames <- NULL
	inDirs <- NULL
	for (d in inDir) {
		message(glue("Searching files with {grepPattern} pattern under {d} ..."))
		tmp <- dir(path=d,pattern=grepPattern)
		message(sprintf("%i files found", length(tmp)))
		fNames <- append(fNames, tmp)
		inDirs <- append(inDirs, rep(d, length(tmp)))
	}
	sampNames <- extract_sampleId(fNames)
	
	### Exclude samples ###
	if (length(excludeSamp)>0){
		message(sprintf("Excluding {%s}",
		                paste(excludeSamp,collapse=",")))
		idx <- which(sampNames %in% excludeSamp)
		sampNames <- sampNames[-idx]
		fNames <- fNames[-idx]
		inDirs <- inDirs[-idx]
	}
	
	### Construct bsseq object ###
	message(sprintf("about to read bismark from \n- %s", 
	                paste(paste(inDirs,fNames,sep="/"), 
	                      collapse = "\n- ")
	                )
	        )
	bs <- bsseq::read.bismark(
		files=paste(inDirs,fNames,sep="/"),
		colData=DataFrame(inDir = inDirs, 
		                  row.names = sampNames
		                  ),
		rmZeroCov=FALSE,
		strandCollapse=TRUE,
		verbose=TRUE 
	)
	
	if (drop_hg38_alt) {
	  # remove spike-in and alt chrs
		bs <- bs[bs@rowRanges@seqnames %in% paste0("chr", c(as.character(1:22), "X", "Y"))] 
	} else {
	  # only remove spike-in
		bs <- bs[bs@rowRanges@seqnames %in% c("pUC19", "J02459.1")] 
	}
	
	### Output ###
	return(list(files=fNames,samples=sampNames,bs=bs))
}


#' Load CpG DMRs of RL-VZ vs SVZ
#' @param inFile (character) The directory of the file
#' @return (GRanges) Return a Granges object of the dmrs
get_cpg_dmrs <- function() {
  inFile <- get_configs("CPG_DMR_FILE")
	message(sprintf("Loading DMRs from: %s", inFile))
	
	df_dmrs <- read.delim(inFile, sep="\t")
	ranges_dmrs <- GRanges(df_dmrs$chr, IRanges(df_dmrs$start, df_dmrs$end))
  mcols(ranges_dmrs)$nCG <- df_dmrs$nCG
  mcols(ranges_dmrs)$diff.Methy <- df_dmrs$diff.Methy
  mcols(ranges_dmrs)$areaStat <- df_dmrs$areaStat

  ### Output ###
	return(ranges_dmrs)
}


#' Load ENCODE cCREs
#' @return (GRanges) Return a Granges object of the ENCODE cCREs
get_encode_cres <- function() {
  inFile <- get_configs("ENCODE_CRE_FILE")
	message(sprintf("Loading ENCODE cCREs from: %s", inFile))

	df_ccres <- read.delim(inFile, header=F, sep="\t")
  
  	ranges_ccres <- GRanges(df_ccres$V1, IRanges(df_ccres$V2, df_ccres$V3))
  	ranges_ccres$ID <- df_ccres$V5
  	ranges_ccres$Type <- df_ccres$V6

  	### Output ###
  	return(ranges_ccres)
}


#' Load RepeatMasker Transposable Elements
#' @return (GRanges) Return a Granges object of the TE
get_repeatmasker_tes <- function() {
  inFile <- get_configs("REPEATMASKER_TE_FILE")
	message(sprintf("Loading RepeatMasker TEs from: %s", inFile))

	ah <- AnnotationHub()
  	query(ah, c("RepeatMasker", "Homo sapiens"))
  	rmskhg38 <- ah[["AH99003"]]

  	### Output ###
  	return(rmskhg38)
}


#' Load GENCODE gene Transcription Start Sites or gene body
#' @param gene_type (character) Either "all" or "protein_coding" genes
#' @param region (character) Either "tss" or "gene_body"
#' @return (GRanges) Return a Granges object of the gene TSS
get_gencode_anno <- function(gene_type = "all", 
                             region = "tss") {
  inFile <- get_configs("GENCODE_GENE_FILE")
	message(glue("Loading GENCODE {gene_type} gene {region} annotation from: {inFile}"))

	### Filter gene type ###
	genes <- rtracklayer::readGFF(inFile)
	if (tolower(gene_type) == "all") {
		genes <- subset(genes, type == "gene")
	} else if (tolower(gene_type) == "protein_coding") {
		genes <- subset(genes, gene_type %in% "protein_coding" & type == "gene")
	} else {
		stop(sprintf("Wrong gene_type provided: %s", gene_type))
	}
	
	### Set ranges ###
	if (tolower(region) == "tss") {
		genes$TSS <- genes$start
		# tss flip for reverse strand
		genes$TSS[which(genes$strand=="-")] <- genes$end[which(genes$strand=="-")] 
		genes$l <- genes$TSS
		genes$r <- genes$TSS
	} else if (tolower(region) == "gene_body") {
		genes$l <- genes$start
		genes$r <- genes$end
	} else {
		stop(sprintf("Wrong region provided: %s", region))
	}

	geneGR <- GRanges(
		genes$seqid,
		IRanges(genes$l, genes$r),
		name=genes$gene_name,
    strand=genes$strand
	) 

	### Output ###
  	return(geneGR)
}


#' infer fetal hindbrain enhancer based on ChIP-seq data
#' @param chipDir (char) Directory of Aldinger ChIP-seq data
#' @param strict (logic) if TRUE, only use intersect of two samples
#' @return (GRanges)
get_fetalCB_enh <- function(chipDir) {
  message(glue("Inferring fetal hindbrain enhancers based on ChIP-seq data ",
               "stored under {chipDir}"))
  
  h3k27ac_1 <- import.bed(sprintf("%s/1_27907-102M_H3K27Ac_hg38.bed",chipDir))
  h3k27ac_2 <- import.bed(sprintf("%s/2_27556-132M_H3K27Ac_hg38.bed",chipDir))
  h3k4me3_1 <- import.bed(sprintf("%s/3_27907-102M_H3K4me3_hg38.bed",chipDir))
  h3k4me3_2 <- import.bed(sprintf("%s/4_27556-132M_H3K4me3_hg38.bed",chipDir))
  
  tss <- get_gencode_anno()
  tss_2kb <- tss
  start(tss_2kb) <- start(tss) - 2000
  end(tss_2kb) <- end(tss) + 2000
  
  ### Combine sample peaks ###
  h3k27ac <- GenomicRanges::intersect(h3k27ac_1, h3k27ac_2)
  h3k4me3 <- GenomicRanges::intersect(h3k4me3_1, h3k4me3_2)
  
  # following definition of ENCODE
  ## Proximal ##
  p_h3k27ac <- h3k27ac[unique(queryHits(findOverlaps(h3k27ac, tss_2kb)))]
  pELS <- GenomicRanges::setdiff(p_h3k27ac, h3k4me3)
  pELS$type <- "feCB_pELS"
  
  ## Distal ##
  dELS <- h3k27ac[-unique(queryHits(findOverlaps(h3k27ac, tss_2kb)))]
  dELS$type <- "feCB_dELS"
  
  ### Final ###
  fetal_CB_enh <- c(pELS, dELS)
  
  return(fetal_CB_enh)
}


get_fetalCB_promoter <- function(chipDir) {
  message(glue("Inferring fetal hindbrain promoters based on ChIP-seq data ",
               "stored under {chipDir}"))

  h3k4me3_1 <- import.bed(sprintf("%s/3_27907-102M_H3K4me3_hg38.bed",chipDir))
  h3k4me3_2 <- import.bed(sprintf("%s/4_27556-132M_H3K4me3_hg38.bed",chipDir))
  
  tss <- get_gencode_anno()
  tss_200 <- tss
  start(tss_200) <- start(tss) - 200
  end(tss_200) <- end(tss) + 200
  
  ### Combine sample peaks ###
  h3k4me3 <- GenomicRanges::intersect(h3k4me3_1, h3k4me3_2)
  
  # following definition of ENCODE
  ### PLS ###
  pls <- h3k4me3[unique(queryHits(findOverlaps(h3k4me3, tss_200)))]
  pls$type <- "feCB_PLS"
  
  return(pls)
}

get_fetalCB_se <- function(dir) {
  message(glue("Collecting fetal hindbrain SE stored under {dir}"))
  
  se <- import.bed(sprintf("%s/superEnhancerLikeElements.bed",dir))
  se$type <- "feCB_SE"
  
  return(se)
}

#' Load Aldinger fetal cerebellar cCREs
#' @return (GRanges) Return a Granges object of the fetal cerebellar cCREs
get_fetalcb_cres <- function() {
  inDir <- get_configs("ALDINGER_CB_CRE_DIR")
	message(glue("Loading Aldinger fetal cerebellar cCREs from: {inDir}"))

	## Load ###
	feCB_enh <- get_fetalCB_enh(glue("{inDir}/raw"))
	feCB_pro <- get_fetalCB_promoter(glue("{inDir}/raw"))
	feCB_se <- get_fetalCB_se(inDir)

  ### Combine cCREs ###
  feCB_CRE <- c(feCB_enh, feCB_pro, feCB_se)
  print(table(feCB_CRE$type))

	### Output ###
  	return(feCB_CRE)
}








