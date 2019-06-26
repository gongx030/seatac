#' vplot
#'

vplot <- function(gr, which, txdb){

  which <- gsub(',', '', which)
  if (!grepl('chr.+:\\d+-\\d+', which))
		stop('which must be in the correct format')
    
  chrom <- gsub('(.+):(.+)-(.+)', '\\1', which)
	start <- as.numeric(gsub('(.+):(.+)-(.+)', '\\2', which))
	end <- as.numeric(gsub('(.+):(.+)-(.+)', '\\3', which))

	which <- GRanges(seqnames = chrom, range = IRanges(start, end))

  if (is.null(mcols(gr)$counts))
    stop('counts field must be available')

  if (is.null(mcols(gr)$predicted_counts))
    stop('predicted_counts field must be available')

  gr <- subsetByOverlaps(gr, which)
  flog.info(sprintf('found %d %d bp bins overlapping with %s:%d-%d', length(gr), metadata(gr)$bin_size, seqnames(which), start(which), end(which)))


	browser()
	tx <- subsetByOverlaps(transcripts(txdb), which)
	gtrack <- GenomeAxisTrack()
	atrack <- AnnotationTrack(gr , name = 'CpG')
	grtrack <- GeneRegionTrack(tx, name = 'Gene Model')
	plotTracks(list(gtrack, atrack, grtrack), from = start, to = end, chromosome = chrom)




	image(mcols(gr)$predicted_counts)
	image(mcols(gr)$counts)

#  plot(mcols(gr)$coverage[, 1])
  image(mcols(gr)$predicted_counts, breaks = c(seq(0, 0.15, length.out = 100), 1), col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
#  image(mcols(gr)$predicted_counts, col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
  y <- summary(as(mcols(gr)$counts, 'dgCMatrix'))
  points(y[, 1] / nrow(mcols(gr)$counts), y[, 2] / ncol(mcols(gr)$counts), pch = 3, cex = 1.25, col = 'black')
  lines(1:length(gr) / length(gr), mcols(gr)$coverage / max(mcols(gr)$coverage), col = 'yellow')
  plot(rowSums(mcols(gr)$predicted_counts[, 1:5]), type = 'l')

#  browser()

#  p.ideo <- Ideogram(gr)
#  p.ideo <- Ideogram(genome = metadata(gr)$providerVersion)
#  p.ideo

} # vplot
