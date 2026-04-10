# show correlation between promoter methylation changes and gene expression
# over multiple upstream distances used for promoter definition

inDir <- "/home/rstudio/isilon/private/projects/FetalHindbrain/EMseq_FETHB2-FETHB3/DMR_plots/251215"

fPfx <- "DMR_vs_DEG_promoters_correlation_summary_minCvg5_"

fList <- dir(inDir, pattern=fPfx)

dat <- list()
ctr <- 1

bonferroniAlpha <- 0.05 / length(fList)
for (f in fList) {
  d <- read.delim(sprintf("%s/%s", inDir, f), header=TRUE, stringsAsFactors=FALSE)
  dat[[ctr]] <- d
  ctr <- ctr + 1
}

dat2 <- do.call(rbind, dat)
dat2 <- subset(dat2, upstream < 6000)

# plot cor_pearson as a function of upstream as a line plot.
# where pval_pearson < bonferroniAlpha, plot an asterix above the point.
# Make font Helvetica, with base size 18
# Add annotation about bonferroniAlpha on the plot
p <- ggplot(dat2, aes(x = upstream/1000, y = cor_pearson)) +
  geom_line(lwd=2) + geom_point(size=4) +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey",lwd=2) +
  theme_minimal(base_size = 24) +
  labs(title = "DMR vs DEG delta w/ changing updist",
       x = "Upstream distance from TSS (kb)",
       y = "Pearson's r") +
  ylim(min(dat2$cor_pearson)-0.05, 0) + xlim(0.8, max(dat2$upstream)/1000)
p <- p + geom_text(
    data = subset(dat2, pval_pearson < bonferroniAlpha),
    aes(label = "*"),
    vjust = -1.5,
    size = 12
    ) + theme(text = element_text(family = "Helvetica"))
p <- p + annotate("text", x = max(dat2$upstream)/1000, y = max(dat2$cor_pearson)+0.1, 
  label = sprintf("* p < %1.2e", bonferroniAlpha), 
  hjust = 1, size=6)
ggsave(p,file=sprintf("%s/DMR_vs_DEG_promoter_methylation_correlation_vs_upstream_distance.pdf", inDir),
  width = 8, height = 6)

  cat(sprintf("Bonferroni corrected alpha for significance: %1.3e\n", bonferroniAlpha))