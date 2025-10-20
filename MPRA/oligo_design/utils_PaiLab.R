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

#' Stitches TFBSs within the given regions by given window size. Only returns 
#' the interval that is the first one with the most n_TFBSs
#' @param x (GRanges) the region to be stitched, with origin col metadata
#' @param target_len (numeric) the target output length of the windows. Each
#' window expands to target_len with given region left boundary as start
#' @param buffer (numeric) minimum left/right distance from TFBS to boundaries
#' @param keep_meta (char) colnames of metadata to preserve 
#' @return (GRanges) stitched TFBS with target_len
stitchSelect_dmr_tfbs <- function(x, 
                                  target_len = 170, buffer = 15,
                                  keep_meta = NULL) {
  x <- sort(x)
  res <- NULL
  max_n <- 0
  
  for (i in 1:length(x)) {
    interval <- x[i]
    start(interval) <- start(interval) - buffer
    # pseudo-end to account for buffer
    end(interval) <- start(interval) + (target_len - buffer) - 1
    ol <- findOverlaps(x, interval, type = "within")
    
    if (length(ol) > max_n) {
      interval$n_TFBS <- length(ol)
      # actual end
      end(interval) <- end(interval) + buffer
      res <- interval
      max_n <- length(ol)
    }
  }
  
  if (! is.null(keep_meta)) {
    for (col in keep_meta) {
      res@elementMetadata[[col]] <- x@elementMetadata[[col]][1]
    }
  }
  
  return(res)
}


#' get targets that overlap subject
#' @param query (GRanges) target gr
#' @param subject (GRanges) subject gr
#' @param subject_name ("char") name of the subject
#' @return (char) targets in "seqnames_start_end" format
get_ol_targets <- function(query, subject, subject_name) {
  targets <- query[queryHits(findOverlaps(query, subject))]$target
  
  message(sprintf(
    "%d targets (duplication kept; %d unique targets) retrieved from %s overlap", 
    length(targets), 
    length(unique(targets)),
    subject_name
  )
  )
  return(targets)
}


#' find the intersect of two genomic regions and preserve the origin from the target_gr
#' @param target_gr (GRanges) the target to be tracked during intersection. If 
#' there is an origin metadata in the target_gr, it will be preserved instead 
#' of generating new origin.
#' @param ol_gr (GRanges)
#' @return (GRanges)
leftIntersect <- function(target_gr, ol_gr) {
  inter <- GenomicRanges::intersect(target_gr, ol_gr)
  overlap <- findOverlaps(target_gr, inter)
  q <- target_gr[queryHits(overlap)]
  s <- inter[subjectHits(overlap)]
  
  if ("origin" %in% colnames(elementMetadata(target_gr))) {
    s$origin <- q$origin
  } else {
    s$origin <- paste(seqnames(q), start(q), end(q), sep = "_")
  }
  
  return(s)
}


#' from a, drop regions having any overlap with b
#' @param a (GRanges) main gr
#' @param b (Granges) overlap gr
#' @param type (char) Definition of overlap; inherited from findOverlaps function.
#' @return (Granges) a without overlapping regions with b
drop_ol <- function(a, b, type = "any") {
  ol <- findOverlaps(a, b, type = type)
  return(a[-queryHits(ol)])
}


#' from a, keep regions having any overlap with b
#' @param a (GRanges) main gr
#' @param b (Granges) overlap gr
#' @param type (char) Definition of overlap; inherited from findOverlaps function.
#' @return (Granges) a with overlapping regions with b
keep_ol <- function(a, b, type = "any") {
  ol <- findOverlaps(a, b, type = type)
  return(a[queryHits(ol)])
}


#' LiftOver from one genome build to another
#' @param gr (GRanges) intervals to liftOver
#' @param chainFile (charaters) path to the chain file
liftOver_gr <- function(gr, 
                        chainFile = glue("/home/rstudio/isilon/src/ucsc-tools",
                                         "/chain_files/hg19ToHg38.over.chain")
) {
  message(sprintf("%i intervals to be lifted.", length(gr)))
  message(sprintf("Chain file: %s", chainFile))
  
  library(liftOver)
  targetBuild <- liftOver(gr, 
                          import.chain(chainFile)
  )
  ln <- unlist(lapply(targetBuild,length))
  if (any(ln!=1)) {
    targetBuild <- targetBuild[-which(ln!=1)] # remove imperfect liftOver mapping
  }
  targetBuild <- unlist(targetBuild)
  
  message(sprintf("%i intervals converted.",length(targetBuild)))
  
  return(targetBuild)
}


#' get promoter regions
#' @param promoter_radius (numeric) Number of base pairs up&down-stream of TSSs. 
#' @return (GRanges)
get_promoters <- function(promoter_radius = 1000) {
  ### set promoter gr ###
  tss <- get_gencode_anno()
  promoter <- tss
  message(sprintf("Promoter is defined as %i bp up- & down- stream of TSS", 
                  promoter_radius
                  )
          )
  start(promoter) <- start(promoter) - promoter_radius
  end(promoter) <- end(promoter) + promoter_radius
  
  return(promoter)
}


#' Remove regions or entire intervals that overlap promoters from provided intervals
#' @param interval (GRanges) Interval to remove promoters
#' @param method (characters) Either "overlap" or "setdiff". Overlap will drop 
#' a entire interval when there is an overlap; while the later one will only remove
#' the region overlapping promoters. Default using overlap. 
#' @param promoter_radius (numeric) Number of base pairs up&down-stream of TSSs. 
#' Default 1kb.
#' @return (GRanges)
drop_promoters <- function(interval, method = "overlap", promoter_radius = 1000) {
  message(sprintf("Received %d intervals", length(interval)))
  ### set promoter gr ###
  promoter <- get_promoters(promoter_radius = promoter_radius)
  
  ### drop ol ###
  if (method == "setdiff") {
    message("Regions overlapping promoters will be dropped.")
    res <- IRanges::setdiff(interval, promoter)
  } else if (method == "overlap") {
    message("Intervals overlapping promoter will be dropped.")
    ol <- findOverlaps(interval, promoter)
    res <- interval[-queryHits(ol)]
  }
  
  message(sprintf("Finally resulted in %i intervals", length(res)))
  return(res)
}


#' get housekeeping genes from msigdb 
#' @param setName (characters) The gene set name under C2 - CGP category
#' @return (character) A vector of housekeeping genes
get_housekeepingGenes <- function(setName = "HOUNKPE_HOUSEKEEPING_GENES") {
  db <- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CGP")
  res <- db %>% 
    filter(gs_name == setName) %>%
    pull(gene_symbol) %>%
    unique()
  message(sprintf("Obtained %i housekeeping genes from %s", 
                  length(res),
                  setName
  )
  )
  
  return(res)
}


#' get neurodevelopment genes
#' @param path (charaters) The path to the file containing the gene list
#' @return (character) A vector of neurodevelopment genes
get_neurodevGenes <- function() {
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
  message(sprintf("Obtained %i neuro development genes from %s", 
                  length(res),
                  paste(c("", paths), collapse = "\n- ")
  )
  )
  
  return(res)
}


#' get G3/4 MB genes
#' @param db (character) The directory to MB_gene list
#' @return (character) A vector of Grp3/4 MB genes
get_g34genes <- function(db = glue("/home/rstudio/isilon/src/gene_list/MB_gene",
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


#' Resize given interval to satisfy certain size requirement
#' @param interval_df (data.frame) A data frame with "seqnames", "start", and 
#' "end" columns.
#' @target_len (numeric) Target length of the new interval. Default 171 bp for MPRA.
#' @param verbose (logical) Output messages.
#' @return (data.frame) with three columns - "seqnames", "start", and "end".
resize_interval <- function(interval, target_len = 171, verbose = T) {
  if (class(interval) == "data.frame") {
    if (! verbose) {message("Using DataFrame mode")}
    width <- interval$end - interval$start + 1 # as IRanges is 1-based
    short <- target_len - width
    left <- floor(short/2)
    right <- short - left
    res <- data.frame(seqnames = interval$seqnames, 
                      start = interval$start - left,
                      end = interval$end + right
    )
  } else if (class(interval) == "GRanges") {
    if (! verbose) {message("Using GRanges mode")}
    width <- width(interval)# as IRanges is 1-based
    short <- target_len - width
    left <- floor(short/2)
    right <- short - left
    res <- interval
    start(res) <- start(interval) - left
    end(res) <- end(interval) + right
  }
  
  return(res)
}


#' Design tiles based on given width and overlap size between tiles
#' Regions will be expanded to exactly fit one or multiple probes 
#' @param gr (GRanges) Input
#' @param width (numeric) Size of target probe
#' @param overlap (numeric) Size of overlap between tiled probes
#' @return (CompressedGRangesList) A list of tiles for each input interval
design_tiles <- function(gr, width = 170, overlap = 50) {
  message(sprintf("Designing tiled probes with width = %d and overlap = %d",
                  width, overlap
                  )
          )
  # <= width expand to width; > width expand to stiched tiles size
  le <- gr[width(gr) <= width]
  gt <- gr[width(gr) > width]
  
  # le
  le <- resize_interval(le, width)
  
  # gt
  n_tiles <- ceiling((width(gt)-width)/(width - overlap)) + 1
  gt <- resize_interval(gt, n_tiles*width - (n_tiles-1)*overlap)
  
  resized_gr <- c(le, gt)
  tiles <- slidingWindows(resized_gr, width = width, step = width - overlap)
  
  # summarize
  stats <- sapply(tiles, length)
  message(
    sprintf("Designed %d tiled probes (median: %d; max: %d) for %d intervals.",
    sum(stats),
    median(stats),
    max(stats),
    length(stats)
    )
    )
  
  return(tiles)
}


#' get ABC predicted CRE targets in Neuronal cells
#' @return (GRanges)
get_neuroABC <- function() {
  res <- import.bed(get_configs("NEURO_ABC_TARGETS"))
  return(res)
}


#' Get Aldinger H3K27ac peak summit
#' @param up (num) length up-stream of summit
#' @param down (num) length down-stream of summit
#' @return (GRanges) regions surrounding summit (May have overlaps)
get_fetalCB_h3k27ac_summit <- function(up = 100, down = 100) {
  message(glue("Selecting for regions {up} bp up- and {down} bp down- stream ",
               "of the peak summit of H3K27ac."
  )
  )
  intervals <- read.csv(glue("/home/rstudio/isilon/private/projects",
                             "/FetalHindbrain/Aldinger_FetalCBL_ChipSeq",
                             "/CBL_Chipseq/FASTQ/cblchipseq/cblchipseq",
                             "/AM SCRI H3K27Ac, H3K4me3, Pol2 MACS2 CSeq 39458",
                             "/017CSCRI_Cerebellum_intervals.csv"), 
                        stringsAsFactors = F, 
                        header = T)
  
  summit <- intervals %>% 
    filter(grepl("H3K27Ac", Sample)) %>%
    mutate(seqnames = sprintf("chr%s", Chromosome),
           start = Peak.Summit - up + 1, 
           end = Peak.Summit + down
    ) %>%
    dplyr::select(seqnames, start, end, Peak.Summit)
  
  summit_gr <- GenomicRanges::makeGRangesFromDataFrame(summit, 
                                                       keep.extra.columns = T
  )
  summit_gr <- unique(summit_gr)
  
  print(quantile(width(summit_gr)))
  
  return(summit_gr)
}


#' check if sequences contain SceI restriction sites for lentiMPRA
#' @param sequences (DNAStringSet) target sequences
check_SceI <- function(sequences) {
  library(Biostrings)
  SceI <- DNAString("TAGGGATAACAGGGTAAT")
  ln <- unlist(lapply(vmatchPattern(SceI, sequences), length))
  
  message(sprintf("%d sequences contain SceI restriction site.", 
                  sum(ln > 0)
                  ))
  
  if (sum(ln > 0)) {
    message(sprintf("Sequences containing SceI: \n%s",
                    paste(names(ln[ln>0]), collapse = "\n")
                    )
            )
  }
}
