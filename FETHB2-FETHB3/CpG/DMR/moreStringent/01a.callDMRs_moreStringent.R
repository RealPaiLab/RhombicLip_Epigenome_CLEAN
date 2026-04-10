# try calling DMRs with more stringent parameters
dmlFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB3/output/downstream/EMseq_FETHB2-FETHB3/CpG/DMRs/CTsnv_excluded/withoutBatchCorrection/240711/DMLs.Rdata"

outDir <- dirname(dmlFile)
outDir <- sprintf("%s/moreStringent_%s", outDir, format(Sys.Date(), "%y%m%d"))
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/callDMRs_moreStringent.log", outDir)
sink(logFile, split=TRUE)
tryCatch({
cat("Loading DMLs from file\n")
t0 <- Sys.time()
load(dmlFile)  # loads dmlTest.sm
print(Sys.time() - t0)

cat("* Calling DMRs with more stringent parameters\n")
params <- list(
    p.threshold = 1e-08,
    minlen = 100,
    minCG = 4,
    dis.merge = 100,
    pct.sig = 0.5,
    delta = 0.1
)

cat(sprintf("Parameters:\n"))
print(params)


cat("Calling DMRs with more stringent parameters\n")
t0 <- Sys.time()
  dmrs <- callDMR(dmlTest.sm,
                  delta=params$delta,
                  p.threshold=params$p.threshold, 
                  minlen=params$minlen,
                  minCG=params$minCG,
                  dis.merge=params$dis.merge,
                  pct.sig=params$pct.sig
  )
print(Sys.time() - t0)
cat(sprintf("Number of DMRs identified: %i\n", nrow(dmrs)))
cat(sprintf("DMR length distribution:"))
print(summary(dmrs$length))

# print histogram plot of DMR diff.Meth, binned by 5%
p <- ggplot(dmrs, aes(x=diff.Methy*100)) +
  geom_histogram(binwidth = 5, color="black", fill="lightblue") +
  labs(title="Histogram of DMR methylation difference",
       x="% DNA methylation difference (RL-VZ - RL-SVZ)",
       y="Number of DMRs") +
  theme_minimal(base_size = 12)
histFile <- sprintf("%s/DMR_diffMethy_histogram_moreStringent.png", outDir)
ggsave(histFile, plot = p, width = 6, height = 4,bg="white")

dmrFile <- sprintf("%s/DMRs_moreStringent.csv", outDir)
write.csv(as.data.frame(dmrs), file = dmrFile, row.names = FALSE)

}, error = function(e) {n
  stop(e)
}, finally = {
  sink()
})