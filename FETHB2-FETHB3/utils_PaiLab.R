library(glue)
library(GenomicRanges)
library(rtracklayer)

#### Utils for DMR analyses ####

#' get preset configs 
#' @param config_name (char) Name of the defined config.
get_configs <- function(config_name = NA) {
  ## CONFIGS ##
  # Cytosine reports #
  CPG_REPORT_DIR <- c(
    glue("/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3",
         "/output/alignment/methyldackel/report/CpG_snpFiltered")
  )
  
  # DMR #
  #CPG_DMR_DATE <- "240522" # outdated
  CPG_DMR_DATE <- "240711"
  CPG_DMR_FILE <- glue(
    "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3",
    "/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded",
    "/withoutBatchCorrection/{CPG_DMR_DATE}/DMRs.csv"
  )
  
  # cCRE #
  ENCODE_CRE_FILE <- glue("/home/rstudio/isilon/private/projects",
                          "/FetalHindbrain/anno/GRCh38-cCREs.bed")
  ALDINGER_CB_CRE_DIR <- glue("/home/rstudio/isilon/private/xsun/output",
                              "/ncMutMB/20240314/cre/Aldinger-FetalCB")
  
  # Gene annotation #
  GENCODE_GENE_FILE <- glue("/home/rstudio/isilon/private/projects/FetalHindbrain",
                            "/anno/gencode.v44.basic.annotation.gtf")
  
  # Transposable elements #
  REPEATMASKER_TE_FILE <- glue("/home/rstudio/isilon/private/xsun/Database",
                               "/RepeatMasker/RepeatMasker_open406_Dec2013_Dfam20_hg38",
                               "/hg38.fa.out.gz")
  
  # metadata #
  METADATA_FILE <- glue("/home/rstudio/isilon/private/projects/FetalHindbrain",
                        "/EMseq_FETHB3/metadata/DNAm_RL_tumours_STables - Table S1.tsv")
  
  # rl meme db #
  RL_MEME_DB <- glue("/data/xsun/db/meme/HOCOMOCOv12/filteredTFdb",
                     "/Hendrikse2022_RL_activeGenes_TF.meme")
  TF_DB <- "HOCOMOCOv12"
  
  # HOCOMOCOv12 
  HOCOMOCO_V12 <- glue("/home/rstudio/isilon/src/transcription-factors",
                       "/HOCOMOCOv12/H12CORE_meme_format.meme")
  
  # Aldinger LCM RL-VZ vs SVZ RNA logFC
  RL_DEG <- glue("/home/rstudio/isilon/private/projects/FetalHindbrain",
                "/AldingerMillen_LCM_RNAseq/output/VZ_SVZ_diffEx/231123",
                "/edgeR_RLVZvsSVZ_231123.txt")
  
  # Neuronal cell ABC targets
  NEURO_ABC_TARGETS <- glue("/data/xsun/20240314/creTargets/Nasser-Neuronal-ABC",
                            "/Nasser-Neuronal-ABC_creTarget_hg38.bed")
  
  # Whalen MPRA results
  WHALEN_MPRA_DIR <- glue("/home/rstudio/isilon/src/neurodev-genomics/MPRA", 
                          "/Whalen_2023/GSE110758_RAW")
  
  # RL active genes
  RL_ACTIVE_GENES <- glue("/data/xsun/db/meme/HOCOMOCOv12/activeGenes",
                          "/Hendrikse2022_RL_activeGenes")
  
  # JASPAR
  JASPAR_DB <- "/home/rstudio/isilon/src/ucsc-goldenpath/hg38/JASPAR2024.bb"
  
  # VISTA positive human hindbrain enhancer
  VISTA_HB <- glue("/home/rstudio/isilon/src/neurodev-genomics",
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
		rmZeroCov=TRUE,
		strandCollapse=TRUE,
		verbose=TRUE 
	)
	
	if (drop_hg38_alt) {
	  # remove spike-in and alt chrs
		bs <- bs[bs@rowRanges@seqnames %in% paste0("chr", c(as.character(1:22), "X", "Y"))] 
	} else {
	  ##### only remove spike-in
		####bs <- bs[bs@rowRanges@seqnames %in% c("pUC19", "J02459.1")] 
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

get_fetalCB_enh <- function() {
    pDir <- "/home/rstudio/isilon/private/xsun/output/ncMutMB/20240314/cre/Aldinger-FetalCB/raw"
    
    ac1 <- import.bed(sprintf("%s/1_27907-102M_H3K27Ac_hg38.bed",pDir))
    ac2 <- import.bed(sprintf("%s/2_27556-132M_H3K27Ac_hg38.bed",pDir))

    ### Combine sample peaks ###
    h3k27ac <- GenomicRanges::intersect(ac1, ac2)

    return(h3k27ac)
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

get_fetalCB_se <- function() {

  inDir <- "/home/rstudio/isilon/private/xsun/output/ncMutMB/20240314/cre/Aldinger-FetalCB"
  
  se <- import.bed(sprintf("%s/superEnhancerLikeElements.bed",inDir))
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

#' get neurodevelopment genes
#' @param path (charaters) The path to the file containing the gene list
#' @return (character) A vector of neurodevelopment genes
get_neurodevGenes <- function(verbose=TRUE) {
  additional_genes = c("OLIG3") # mentioned by Kim Aldinger
  
  paths <- c(avc = glue("/home/rstudio/isilon/src/gene_list/brain-development",
                        "/aldinger2021_vladoiu2019_carter2018.csv"),
             leto = glue("/home/rstudio/isilon/src/gene_list/brain-development",
                  "/LetoEtAl.txt"),
             bhaduri = glue("/home/rstudio/isilon/src/gene_list/brain-development",
                            "/BhaduriKriegstein2021_nature_neocortex.txt"),
             sepp = glue("/home/rstudio/isilon/src/gene_list/brain-development",
                         "/Sepp2023_CB.txt"),
             ian_leo = glue("/home/rstudio/isilon/src/gene_list/brain-development",
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
  if (verbose) {
    message(sprintf("Obtained %i neuro development genes from %s", 
                    length(res),
                    paste(c("", paths), collapse = "\n- ")
    )
    )
  }
  
  return(res)
}

#' get G3/4 MB genes
#' @param db (character) The directory to MB_gene list
#' @return (character) A vector of Grp3/4 MB genes
get_g34genes <- function(db = glue("/home/rstudio/isilon/src/gene_list/MB_gene",
                                   "/MBgene_database_20240709171001.csv"),
                                   verbose = TRUE                            
                         ) {
  df <- read.csv(db, stringsAsFactors = F, header = T, row.names = 1)
  sub <- df[,grepl("Hendrikse2022_G3_G4_MB_genes|Northcott2017_G34_genes", 
                   colnames(df)
                   )
            ]

  g34_genes <- unique(rownames(sub[rowSums(sub) != 0,]))
  
  if (verbose) {
    message(sprintf("Obtained %i G3/4 MB genes from %s", 
                    length(g34_genes),
                    db
    ))
  }
  return(g34_genes)
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


#' get nearest gene to all ranges in gr
getNearestGene <- function(gr, gene_types = NULL,verbose=FALSE){
    geneFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v42.basic.annotation.gtf"
    
    if (verbose) cat("Reading gene annotation...\n")
    genes <- rtracklayer::readGFF(geneFile)
    genes <- subset(genes, type == "gene")

    if (! is.null(gene_types)) {
      genes <- subset(genes, gene_type %in% gene_types)
    }
    if (verbose) cat(sprintf("Got %i genes of %s type\n", nrow(genes), (gene_types)))

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

#' get targets that overlap subject
#' @param query (GRanges) target gr
#' @param subject (GRanges) subject gr
#' @param subject_name ("char") name of the subject
#' @return (char) targets in "seqnames_start_end" format
get_ol_targets <- function(query, subject, subject_name) {
  targets <- query[queryHits(findOverlaps(query, subject))]
  
  message(sprintf(
    "%d targets (duplication kept; %d unique targets) retrieved from %s overlap", 
    length(targets), 
    length(unique(targets)),
    subject_name
  )
  )
  return(targets)
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

