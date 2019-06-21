#' vplot
#'

setGeneric('vplot', function(gr, which, ...) standardGeneric('vplot'))

setMethod(
  'vplot',
  signature = c(
    gr = 'GRanges',
    which = 'character'
  ),
  function(gr, which, ...){

    which <- gsub(',', '', which)
    if (!grepl('chr.+:\\d+-\\d+', which))
      stop('which must be in the correct format')
    
    chrom <- gsub('(.+):(.+)-(.+)', '\\1', which)
    start <- as.numeric(gsub('(.+):(.+)-(.+)', '\\2', which))
    end <- as.numeric(gsub('(.+):(.+)-(.+)', '\\3', which))

    which <- GRanges(seqnames = chrom, range = IRanges(start, end))

    vplot_core(gr, which)
  }
)


vplot_core <- function(gr, which){

  if (is.null(mcols(gr)$counts))
    stop('counts field must be available')

  if (is.null(mcols(gr)$predicted_counts))
    stop('predicted_counts field must be available')

  gr <- subsetByOverlaps(gr, which)
  flog.info(sprintf('found %d %d bp bins overlapping with %s:%d-%d', length(gr), metadata(gr)$bin_size, seqnames(which), start(which), end(which)))

  browser()
#  plot(mcols(gr)$coverage[, 1])
  image(mcols(gr)$predicted_counts, col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
  y <- summary(as(mcols(gr)$counts, 'dgCMatrix'))
  points(y[, 1] / nrow(mcols(gr)$counts), y[, 2] / ncol(mcols(gr)$counts), pch = 3, cex = 1.25, col = 'black')
#  lines(1:length(gr) / length(gr), mcols(gr)$coverage[, 1] / max(mcols(gr)$coverage[, 1]), col = 'yellow')
  lines(1:length(gr) / length(gr), mcols(gr)$coverage[, 2] / max(mcols(gr)$coverage[, 2]), col = 'yellow')

  p.ideo <- Ideogram(gr)
  p.ideo <- Ideogram(genome = metadata(gr)$providerVersion)
  p.ideo

  browser()
} # vplot
