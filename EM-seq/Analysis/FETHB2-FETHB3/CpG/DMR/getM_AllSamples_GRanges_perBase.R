# computes methylation distribution across GenomicRanges. Returns per-base measure.
source("getRec_GRanges.R")
suppressMessages(require(dplyr))
 
#' @param (char) cytoDir path to cytosine report files
#' @param (GRanges) loci for which %M is desired
#' @return (list) 2 tables, rows are loci, columns are samples
#' 1. pctM: % methylation
#' 2. cov: coverage
getM_AllSamples_GRanges_perBase <- function(cytoDir, gr, incSamples=NULL){
fList <- dir(path=cytoDir,pattern="cytosine_report.txt.gz$")

if (!is.null(incSamples)){
	sampName <- sub(".cytosine_report.txt.gz","",fList)
	fList <- fList[which(sampName %in% incSamples)]
}
cat(sprintf("Got %i samples\n",length(fList)))

out <- list()
for (samp in fList) {  #loop over samples
    sampName <- sub(".cytosine_report.txt.gz","",samp)
	inFile <- sprintf("%s/%s",cytoDir,samp)
	cat("---------------------------\n")
	cat(sprintf("%s\n",samp))
	cat("---------------------------\n")
	baseF <- sub("cytosine_report.txt.gz","",basename(samp))
	
	# get methylation counts related to each target
	rec <- getRec_GRanges(inFile,gr,verbose=FALSE,numCores=16)
	y <- lapply(rec, function(x) {
		if (length(dim(x))<2) return(cbind(C=0,COV=0,pctM=NA))
		else {
			x$pos <- as.integer(x$pos)
			x$pos[which(x$strand == "-")] <- x$pos[which(x$strand == "-")]-1
			y <- x %>% group_by(chr,pos) %>%
				summarise(C_count=sum(C_count), T_count=sum(T_count))
			
			origx <- x
			x <- y
            cov <- x$C_count + x$T_count
			z <- cbind(
                    C=x$C_count,
                    COV=cov,
				    pctM=x$C_count/cov
			)
            z
            }
    })
			#})
	names(y) <- gr$name
	
    out[[sampName]] <- y
	rm(y)
}

return(out)
}