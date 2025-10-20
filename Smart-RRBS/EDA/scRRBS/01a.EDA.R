library(ggplot2)
library(methylKit)


rrbs_root = "/.mounts/labs/pailab/private/projects/MB_scDNAme/input/SmartRRBS_pilot/raw/RRBS" # Sites matrices main folder
matrix_files <- list.files(rrbs_root, recursive = T, pattern = "site_matrix.csv")
print(matrix_files)

### matrix size checked by command line ###
# wc -l ./*/site_matrix.csv
# 10374318 ./1170_1/site_matrix.csv
# 2486428 ./1170_2_CD/site_matrix.csv
# 9853496 ./1170_2/site_matrix.csv

# ls ./*/site_matrix.csv | while read line; do echo -ne $line"\t"; head -n1 $line | awk 'BEGIN{FS=","} {print NF}'; done
# ./1170_1/site_matrix.csv	97
# ./1170_2_CD/site_matrix.csv	24
# ./1170_2/site_matrix.csv	73

### Check one matrix ###
if (F) { # not needed anymore
  ## use the smallest matrix as input for exploration ##
  test_file <- sprintf("%s/%s", rrbs_root, "1170_2_CD/site_matrix.csv")
  test_mat <- read.csv(test_file, header = T, stringsAsFactors = F, na.strings = "-", row.names = 1)
  print(test_mat[1:3,1:3])
  
  ## Check sample and site level coverage ##
  barplot(sort(colSums(!is.na(test_mat))), las = 2)
  perSample_cov <- colSums(!is.na(test_mat))/nrow(test_mat)
  
  print(quantile(colSums(!is.na(test_mat))))
  ggplot(as.data.frame(perSample_cov), aes(x = "", y = perSample_cov)) + geom_violin() + geom_jitter() + theme_minimal()
  
  
  ### load with MethylKit and merge regions
  cov_dir <- sprintf("%s/%s", rrbs_root, "1170_2_CD/cov_files")
  
  myobj <- methRead(as.list(list.files(cov_dir, pattern = ".bismark.cov.gz", full.names = T)),
                    sample.id = as.list(list.files(cov_dir, pattern = ".bismark.cov.gz", full.names = F)),
                    assembly = "hg38",
                    treatment = rep(0, length(list.files(cov_dir, pattern = ".bismark.cov.gz", full.names = T))),
                    context = "CpG",
                    pipeline = "bismarkCoverage",
                    mincov = 1,
                    header = F
  )
  tiles <- tileMethylCounts(myobj, win.size=100000,step.size=100000,cov.bases = 5, mc.cores = 20)
  print(sort(unlist(lapply(tiles, nrow))))
  
  keptCells <- unlist(lapply(tiles, function(x) {if (nrow(x) > 100) {return(x@sample.id)}}))
  
  filtered_tiles <- reorganize(tiles,sample.ids=keptCells, treatment = rep(0, length(keptCells)))
  message(sprintf("%d Cells kept after removing low coverage cells", length(keptCells)))
  
  meth <- unite(filtered_tiles, destrand=T, min.per.group = as.integer(length(filtered_tiles)*1))
  tmp <- percMethylation(meth)
  print(dim(tmp))
  
  heatmap(tmp)
  backup <- tiles
}


### Load all three plates together with MethylKit
# the overlapping regions are just too low for this plate. Having the other two better quality plates may give better results

cov_dir <- paste(rrbs_root, c("1170_1", "1170_2_CD", "1170_2"), "cov_files", sep = "/")
cells <- list.files(cov_dir, pattern = ".bismark.cov.gz", full.names = F)
myobj <- methRead(as.list(list.files(cov_dir, pattern = ".bismark.cov.gz", full.names = T)),
                  sample.id = as.list(cells),
                  assembly = "hg38",
                  treatment = ifelse(startsWith(cells, "1170_1"), 0, 1),
                  context = "CpG",
                  pipeline = "bismarkCoverage",
                  mincov = 1,
                  header = F
)

## merge input windows
tiles <- tileMethylCounts(myobj, win.size=10000,step.size=10000,cov.bases = 1, mc.cores = 20)
n_windows <- unlist(lapply(tiles, nrow))
print(sort(n_windows))

hist(n_windows, 100)
boxplot(n_windows ~ unlist(lapply(tiles, function(x){stringr::str_extract_all(x@sample.id, pattern = "1170_\\d(_CD)?")})), xlab = "plate")

## remove cells with low window coverage
quant <- 0.6
message(sprintf("%f quantile window num: %f", quant, quantile(n_windows, quant)))
window_num_threshold <- quantile(n_windows, quant)
keptCells <- unlist(lapply(tiles, function(x) {if (nrow(x) >= window_num_threshold) {return(x@sample.id)}}))

filtered_tiles <- reorganize(tiles,sample.ids=keptCells, treatment = rep(0, length(keptCells)))
message(sprintf("%d Cells kept after removing low coverage cells from a total of %d cells", length(keptCells), length(cells)))
n_windows <- unlist(lapply(filtered_tiles, nrow))
print(sort(n_windows))
boxplot(n_windows ~ unlist(lapply(filtered_tiles, function(x){stringr::str_extract_all(x@sample.id, pattern = "1170_\\d(_CD)?")})), xlab = "plate")

## combine all cells and keep windows covered by all cells
meth <- unite(filtered_tiles, min.per.group = as.integer(min(table(filtered_tiles@treatment))*1), mc.cores = 20)
tmp <- percMethylation(meth)
colnames(tmp) <- stringr::str_split(colnames(tmp), pattern = "\\.R1_bismark_bt2_pe", simplify = T)[,1]
rownames(tmp) <- paste0(meth$chr, ":", meth$start, "-", meth$end)
print(dim(tmp))

## remove low variance windows? 
row_variances <- apply(tmp, 1, function(x) {var(x, na.rm = T)})
high_var <- names(row_variances[row_variances >= quantile(row_variances, 0.)]) # used all windows as the number was low
print(length(high_var))

## plot heatmap
anno <- data.frame(plate = unlist(stringr::str_extract_all(colnames(tmp), pattern = "1170_\\d(_CD)?")))
rownames(anno) <- colnames(tmp)
pheatmap::pheatmap(tmp[high_var,], cluster_rows = T, cluster_cols = T, 
                   na_col = "green", annotation_col = anno, border_color = NA,
                   show_colnames = F, show_rownames = T)

## plot PCA
pcs <- PCASamples(meth, obj.return = T)
pcs$x
df <- as.data.frame(pcs$x)
df$plate <- unlist(stringr::str_extract_all(rownames(df), pattern = "1170_\\d(_CD)?"))
ggplot(df, aes(x = PC1, y = PC2, color = plate)) + geom_point(size = 3) + theme_minimal()



















