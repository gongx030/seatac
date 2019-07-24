#' readFragmentSizeMatrix
#' @export
#'
readFragmentSizeMatrix <- function(
	x, 
	windows, 
	window_size = 500,
	bin_size = 10,
	fragment_size_range = c(50, 670), 
	fragment_size_interval = 20
){

  num_samples <- metadata(x)$num_samples
  n_bins_per_window <- window_size / bin_size
	breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)

	n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval + 1

	windows <- windows %>% resize(fix = 'center', width = window_size)

	wb <- cbind(rep(1:length(windows), n_bins_per_window), 1:(length(windows)* n_bins_per_window))	# windows ~ bins, sorted by window

	bins <- windows %>%
		slidingWindows(width = bin_size, step = bin_size) %>%
		unlist()

	# coverage of each 10-bp bin in [-expand / 2, +expand / 2] window
	gr <- rep(windows, num_samples)
	mcols(gr)$group <- rep(1:num_samples, each = length(windows))

	mcols(gr)$coverage <- do.call('rbind', lapply(1:num_samples, function(i){
		cvg <- coverage(x[mcols(x)$group == i])
		matrix(mean(cvg[bins]), nrow = length(windows), ncol = n_bins_per_window, byrow = TRUE)
	}))

	# compute the center point between PE reads
	# this is faster than using GAlignmentPairs
	x <- x[strand(x) == '+']
	x <- GRanges(
		seqnames = seqnames(x), 
		range = IRanges(start(x) + round(mcols(x)$isize / 2), width = 1), 
		isize = mcols(x)$isize,
		group = mcols(x)$group
	)

  x$fragment_size <- as.numeric(cut(x$isize, breaks))	# discretize the fragment size
  x <- x[!is.na(x$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"

  X <- Reduce('rbind', lapply(1:num_samples, function(i){

	  # find the # PE read centers per window per sample
	  xi <- subsetByOverlaps(x[mcols(x)$group == i], windows)

    CF <- sparseMatrix(i = 1:length(xi), j = xi$fragment_size, dims = c(length(xi), length(breaks)))  # read center ~ fragment size
    BC <- as.matrix(findOverlaps(bins, xi))	# bins ~ read center
    BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(xi))), 'dgCMatrix') # bins ~ read center
    BF <- BC %*% CF  # bins ~ fragment size
    BF[BF > 0] <- 1
    BF <- as.matrix(BF[wb[, 2], ])
    dim(BF) <- c(n_bins_per_window, length(windows), n_intervals)	# convert BF into an array with bathc_size ~ n_bins_per_window ~ n_intervals
    BF <- aperm(BF, c(1, 3, 2)) # n_bins_per_window ~ n_intervals ~ bathc_size 
    dim(BF) <- c(n_bins_per_window * n_intervals, length(windows))
    BF <- aperm(BF, c(2, 1))  # bathc_size ~ n_bins_per_window * n_intervals
    BF <- as(BF, 'dgCMatrix')
    BF
  }))

  mcols(gr)$counts <- X

	mcols(gr)$min_coverage <- rowMin(mcols(gr)$coverage)
	mcols(gr)$max_coverage <- rowMax(mcols(gr)$coverage)
	mcols(gr)$coverage <- (mcols(gr)$coverage - mcols(gr)$min_coverage) / (mcols(gr)$max_coverage - mcols(gr)$min_coverage)
	mcols(gr)$num_reads <- rowSums(X)

  metadata(gr)$fragment_size_range  <- fragment_size_range
  metadata(gr)$fragment_size_interval <- fragment_size_interval
  metadata(gr)$bin_size <- bin_size
  metadata(gr)$window_size <- window_size 
  metadata(gr)$n_intervals <- n_intervals
  metadata(gr)$num_samples <- num_samples
  metadata(gr)$n_bins_per_window <- n_bins_per_window 


	gr

} # readFragmentSizeMatrix
