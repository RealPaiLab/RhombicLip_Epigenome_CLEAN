library(Gviz)

# Sample data for the custom track
data <- data.frame(
  start = c(10, 50, 100, 150),
  end = c(20, 60, 110, 160),
  score = c(5, 15, 10, 20)
)

# Custom plotting function
plottingFunction <- function(GdObject, prepare = FALSE, ...) {
  if (prepare) {
    return(GdObject)
  } else {
    dp <- displayPars(GdObject)    
    dat <- GdObject@variables
    plot(dat$start, dat$score,
         type = "h",
         col = dp$col,
         lwd = dp$lwd,
         xlab = "Position",
         ylab = "Score",
         main = "blah"
    )
  }
  return(GdObject)
}

# Create the CustomTrack object
track <- CustomTrack(
  plottingFunction = plottingFunction,
  variables = data,
  name = "My Custom Track",
  genome = "hg19",
  chromosome = "chr1",
  start = 0,
  end = 200,
  col = "red",
  lwd = 2
)

# Plot the track
pdf("gviz_test.pdf")
plotTracks(track)
dev.off()