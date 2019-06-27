#' vplot
#'

vplot <- function(gr, which, txdb){


	if (is.character(which)){
	  which <- gsub(',', '', which)
 	 	if (!grepl('chr.+:\\d+-\\d+', which))
			stop('which must be in the correct format')
    
		chrom <- gsub('(.+):(.+)-(.+)', '\\1', which)
		start <- as.numeric(gsub('(.+):(.+)-(.+)', '\\2', which))
		end <- as.numeric(gsub('(.+):(.+)-(.+)', '\\3', which))

		which <- GRanges(seqnames = chrom, range = IRanges(start, end))
	}

  if (is.null(mcols(gr)$counts))
    stop('counts field must be available')

  if (is.null(mcols(gr)$predicted_counts))
    stop('predicted_counts field must be available')

	bin_size <- 10
	start(which) <- floor(start(which) / bin_size) * bin_size + 1
	end(which) <- ceiling(end(which) / bin_size) * bin_size 
  bins <- slidingWindows(which, width = bin_size, step = bin_size)
	bins <- Reduce('c', bins)

  gr <- subsetByOverlaps(gr, which)

	A <- as.matrix(findOverlaps(bins, gr))
	A <- sparseMatrix(i = A[, 1], j = A[, 2], dims = c(length(bins), length(gr)))
	X <- A %*% mcols(gr)$counts %>% as.matrix()
	Xp <- A %*% mcols(gr)$predicted_counts %>% as.matrix()
	valid <- rowSums(A) > 0
	X[!valid, ] <- NA
	Xp[!valid, ] <- NA
	C <- A %*% mcols(gr)$cluster %>% as.matrix()

	genome_coords <- format(seq(start(which), end(which), length.out = 5), nsmall = 0, big.mark = ',')
	coords <- seq(0, 1, length.out = 5)

	fs_coords2 <- seq(100, 600, by = 100)

	par(mfrow = c(3, 1))

	par(mar = c(2, 5, 2, 5))
	image(X, axes = FALSE, col = c('blue', 'red'), xaxt = 'n')
	axis(side = 1, at = coords, labels = genome_coords, tick = TRUE, cex.axis = 1.5)
	axis(side = 2, at = (fs_coords2 - 50) / (650 - 50), label = fs_coords2, tick = TRUE, cex.axis = 1.5, las = 2)

	par(mar = c(2, 5, 1, 5))
	image(Xp, axes = FALSE, breaks = c(0, seq(0, quantile(Xp, 0.95, na.rm = TRUE), length.out = 99), 1), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'))
	axis(side = 1, at = coords, labels = genome_coords, tick = TRUE, cex.axis = 1.5)
	axis(side = 2, at = (fs_coords2 - 50) / (650 - 50), label = fs_coords2, tick = TRUE, cex.axis = 1.5, las = 2) 

	nfr <- rowMeans(Xp[, 1:((110 - 50) / 20)])
	mono <- rowMeans(Xp[, c(round((180 - 50) / 20):round((247 - 50) / 20))])

	par(mar = c(3, 5, 1, 5))
	par(xaxs = 'i') 
	xp <- seq(0, 1, length.out = nrow(X))
	plot(xp, nfr, type = 'l', col = 'black', lwd = 2, axes = FALSE, ylim = c(0, 0.4), ylab = '')
	lines(xp, mono, type = 'l', col = 'red', lwd = 2)
	axis(side = 1, at = coords, labels = genome_coords, tick = TRUE, cex.axis = 1.5)
	axis(side = 2, at = seq(0, 0.4, by = 0.1), tick = TRUE, cex.axis = 1.5, las = 2)
	legend('topright', c('NFR', 'Mono-nucleosome'), lwd = 2, col = c('black', 'red'), cex = 1.5, bty = 'n')

	browser()



	tx <- subsetByOverlaps(transcripts(txdb), which)

	di = rowMeans(Xp[, round((315 - 50) / 20):round((473 - 50) / 20)])
	tri = rowMeans(Xp[, round((558 - 50) / 20):round((615 - 50) / 20)])



#	image(C, axes = FALSE, col = c('white', 'black'))
#	ch <- 1 / (ncol(C) - 1)
#	cw <- 1 / (nrow(C) - 1)
#	hline <- seq(-ch / 2, 1 + ch / 2, length.out = ncol(C) + 1)
#	abline(h = hline, col = 'black', lwd = 1, xpd = FALSE)	
#	abline(v = c(0, 1), col = 'black', lwd = 1)
	

	browser()


	p <- autoplot(txdb, which = which)
	se <- SummarizedExperiment(rowRanges = bins, assays = list(counts = X))
	p2 <- autoplot(se, geom = 'line')
	tracks(p, p2)

	gtrack <- tracks('Gene' = p, xlim = which)
	gtrack <- as(gtrack, 'grob')
	counts_track <- autoplot(t(X))
	grid.arrange(gtrack,counts_track, predicted_counts_track)

	predicted_counts_track <- autoplot(t(Xp))
	grid.arrange(gtrack,counts_track, predicted_counts_track)

	p2 <- tracks(Counts = p2, xlim = which)
	p2 <- as(p2, 'grob')
	grid.arrange(gtrack, p2)


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

add.grid <- function(x, col = 'gray', lwd = 0.1, v = TRUE, h = TRUE){

	if (!is.matrix(x))
		stop('x must be a matrix')

	if (ncol(x) == 1)
		ch <- 1
	else
		ch <- 1 / (ncol(x) - 1)

	if (nrow(x) > 1)
		cw <- 1 / (nrow(x) - 1)
	
	vline <- NULL
	hline <- NULL

	if (v){
		if (nrow(x) == 1)
			vline <- c(-1, 1)
		else
			vline <- seq(-cw / 2, 1 + cw / 2, length.out = nrow(x) + 1)
		abline(v = vline, col = col, lwd = lwd, xpd = FALSE)
	}

	if (h){
		if (ncol(x) == 1)
			hline <- c(-1, 1)
		else
			hline <- seq(-ch / 2, 1 + ch / 2, length.out = ncol(x) + 1)
		abline(h = hline, col = col, lwd= lwd, xpd = FALSE)	
	}
	list(vline = vline, hline = hline)	
}

