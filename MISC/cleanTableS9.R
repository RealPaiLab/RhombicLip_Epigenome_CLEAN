# remove Mannens ref and GISTIC OL columns from Table S9 as we have better sources now.
rm(list=ls())

library(readxl)


suppFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/ms_SuppTable/RhombicLipEpigenome_SuppTables_260331.xlsx"
outDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/fixTableS9"
if (!file.exists(outDir)) dir.create(outDir, recursive=FALSE)

 cat("Starting DMR CNA overlap analysis...\n")
df <- as.data.frame(read_excel(suppFile, sheet="TableS9"))

idx <- c(which(colnames(df) == "Mannens_gene"),grep("GISTIC",colnames(df)))
cat(sprintf("Removing columns: %s\n", paste(colnames(df)[idx], collapse=", ")))
df <- df[, -idx]

cat("removing duplicate rows\n")
cat(sprintf("before: %i rows\n", nrow(df)))
df <- df[!duplicated(df),]
cat(sprintf("after: %i rows\n", nrow(df)))

write.table(df, file=sprintf("%s/TableS9_clean.txt", outDir), sep="\t", quote=FALSE, row.names=FALSE)
