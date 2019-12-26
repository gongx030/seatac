#' read_fragment_size_matrix
#'
#' @param x a GAlignments object.
#' @param windows a GRange object that define a set of genomic regions.
#' @param window_size The size of each genomic window for training the model.
#' @param bin_size The bin size
#' @param fragment_size_range fragment_size_range
#' @param fragment_size_interval fragment_size_interval
#'
#' @export
#'
read_fragment_size_matrix <- function(
	x, 
	windows, 
	window_size = 320,
	bin_size = 5,
	fragment_size_range = c(50, 690), 
	fragment_size_interval = 10
){

	windows <- resize(windows, fix = 'center', width = window_size)

  num_samples <- metadata(x)$num_samples

  n_bins_per_window <- window_size / bin_size

	breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
	centers <- (breaks[-1] + breaks[-length(breaks)]) / 2

	n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval

	flog.info(sprintf('fragment size range(fragment_size_range):%d-%d', fragment_size_range[1], fragment_size_range[2]))
	flog.info(sprintf('fragment size interval(fragment_size_interval):%d', fragment_size_interval))
	flog.info(sprintf('# intervals:%d', n_intervals))

	wb <- cbind(rep(1:length(windows), n_bins_per_window), 1:(length(windows)* n_bins_per_window))	# windows ~ bins, sorted by window

	bins <- windows %>%
		slidingWindows(width = bin_size, step = bin_size) %>%
		unlist()

	gr <- rep(windows, num_samples)
	mcols(gr)$group <- rep(1:num_samples, each = length(windows))

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

    CF <- sparseMatrix(i = 1:length(xi), j = xi$fragment_size, dims = c(length(xi), n_intervals))  # read center ~ fragment size
    BC <- as.matrix(findOverlaps(bins, xi))	# bins ~ read center
    BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(xi))), 'dgCMatrix') # bins ~ read center
    BF <- BC %*% CF  # bins ~ fragment size
    BF <- as.matrix(BF[wb[, 2], ])
    dim(BF) <- c(n_bins_per_window, length(windows), n_intervals)	# convert BF into an array with n_bins_per_window ~ batch_size ~ n_intervals
    BF <- aperm(BF, c(2, 1, 3)) # batch_size, n_bins_per_window ~ n_intervals
    dim(BF) <- c(length(windows), n_bins_per_window * n_intervals)
    BF <- as(BF, 'dgCMatrix')
    BF
  }))

  mcols(gr)$counts <- X
	mcols(gr)$num_reads <- Matrix::rowSums(X)
	mcols(gr)$window_id <- rep(1:length(windows), num_samples)

	metadata(gr)$nfr <- centers <= 100
	metadata(gr)$mono_nucleosome <- centers >= 180 & centers <= 247

  metadata(gr)$fragment_size_range  <- fragment_size_range
  metadata(gr)$fragment_size_interval <- fragment_size_interval
  metadata(gr)$bin_size <- bin_size
  metadata(gr)$window_size <- window_size 
  metadata(gr)$n_intervals <- n_intervals
  metadata(gr)$num_samples <- num_samples
  metadata(gr)$n_bins_per_window <- n_bins_per_window 
  metadata(gr)$breaks <- breaks
  metadata(gr)$centers <- centers

	gr

} # read_fragment_size_matrix

