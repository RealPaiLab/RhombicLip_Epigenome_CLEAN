# computes methylation distribution across GenomicRanges given
source("getRec_GRanges.R")

#' Given R GRanges and S tables, return S by R table of locus-level 
#' percent M for all samples. 
#' @param (char) cytoDir path to cytosine report files
#' @param (GRanges) loci for which %M is desired
#' @return (list) 2 tables, rows are loci, columns are samples
#' 1. pctM: % methylation
#' 2. cov: coverage
getM_AllSamples_GRanges <- function(cytoDir, gr, incSamples=NULL){
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
            cov <- sum(x$C_count) + sum(x$T_count)	
			return(cbind(
                    C=sum(x$C_count),
                    COV=cov,
				    pctM=sum(x$C_count)/cov
            ))}
    })
			#})
	names(y) <- gr$name
	z <- data.frame(do.call("rbind",y))
	z$pctM <- round(z$pctM,digits=3)
	z$target_length <- width(gr)
	z$name <- names(y)

    out[[sampName]] <- z
	rm(y)
}

cvg <- matrix(NA,nrow=length(out),ncol=length(gr))
rownames(cvg) <- names(out);
colnames(cvg) <- names(gr)
pctM <- matrix(NA,nrow=length(out),ncol=length(gr))
rownames(pctM) <- names(out);
colnames(pctM) <- names(gr)


for (k in 1:length(out)){
	tmp <- out[[k]]
	cvg[k,] <- t(tmp$COV)
	pctM[k,] <- t(tmp$pctM)
	
}

return(list(COV=cvg,pctM=pctM))
}