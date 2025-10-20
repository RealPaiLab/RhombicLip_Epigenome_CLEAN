


#' @param gr (GRanges) ranges to check
#' @returns (matrix) rows are regions, columns are TFs, cells count num TFBS
getTFBSMotifMatrix <- function(gr, jaspar_hg38, tmpDir=".", outFile="jaspar.out", convertMat=F){
    grbed <- sprintf("%s/gr.bed", tmpDir)
    export.bed(gr,grbed)
    
    extractTFBS <- "./extract_TFBSs_JASPAR.sh"

    if (file.exists("jaspar.out")) unlink("jaspar.out")
##
    args = c(sprintf("-i %s", grbed),
        sprintf("-b %s", jaspar_hg38),
        sprintf("-o %s -s 400 -p 12", outFile)
    )
    system2(extractTFBS, args)

    if (convertMat) {
        message("converting to matrix")
        mat <- jaspar2matrix(gr,outFile)
        return(mat)
    }
    message(sprintf("--- Finsihed. TFBS outFile stored under %s", outFile))
}

#' convert output of extract_TFBSs_JASPAR to a matrix
#' rows are genomic ranges, columns are TF names. cell values indicate
#' number of TFBS occurrences in that interval.
#' @param gr (GRanges) intervals of interest
#' @param jaspar (char) path to output bed file of extract_TFBSs_JASPAR.sh
jaspar2matrix <- function(gr,jaspar){
    tf <- read.delim(jaspar,h=F)
tfgr <- GRanges(tf[,1],IRanges(tf[,2],tf[,3]))

TFs <- unique(tf[,7])
ol <- findOverlaps(gr,tfgr)
mat <- matrix(0,nrow=length(gr), ncol=length(TFs))
colnames(mat) <- TFs
for (k in unique(queryHits(ol))){
    cur <- subjectHits(ol)[which(queryHits(ol) ==k)]
    cur <- table(tf[cur,7])
    cts <- as.integer(cur)
    for (nm in names(cur)) {
        mat[k,which(colnames(mat)==nm)] <- cts[which(names(cur) == nm)]
    }
}
    return(mat)
}