#' readFragmentSizeMatrix
#'
#' @param x a GAlignments object.
#' @param windows a GRange object that define a set of genomic regions.
#' @param window_size The size of each genomic window for training the model.
#' @param step_size The step size for the moving windows
#' @param bin_size The bin size
#' @param fragment_size_range
#' @param fragment_size_interval
#'
#' @export
#'
readFragmentSizeMatrix <- function(
	x, 
	windows, 
	window_size = 320,
	step_size = 32,
	bin_size = 5,
#	fragment_size_range = c(50, 690), 
#	fragment_size_interval = 10
	fragment_size_range = c(60, 700), 
	fragment_size_interval = 20
){

  num_samples <- metadata(x)$num_samples

	expand <- 2 * window_size - step_size

  n_bins_per_window <- expand / bin_size

	breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)

	n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval

	windows <- windows %>% resize(fix = 'center', width = expand)

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

    CF <- sparseMatrix(i = 1:length(xi), j = xi$fragment_size, dims = c(length(xi), n_intervals))  # read center ~ fragment size
    BC <- as.matrix(findOverlaps(bins, xi))	# bins ~ read center
    BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(xi))), 'dgCMatrix') # bins ~ read center
    BF <- BC %*% CF  # bins ~ fragment size
    BF[BF > 0] <- 1
    BF <- as.matrix(BF[wb[, 2], ])
    dim(BF) <- c(n_bins_per_window, length(windows), n_intervals)	# convert BF into an array with n_bins_per_window ~ batch_size ~ n_intervals
    BF <- aperm(BF, c(2, 1, 3)) # batch_size, n_bins_per_window ~ n_intervals
    dim(BF) <- c(length(windows), n_bins_per_window * n_intervals)
    BF <- as(BF, 'dgCMatrix')
    BF
  }))

  mcols(gr)$counts <- X
	mcols(gr)$num_reads <- rowSums(X)
	mcols(gr)$window_id <- rep(1:length(windows), num_samples)

  metadata(gr)$fragment_size_range  <- fragment_size_range
  metadata(gr)$fragment_size_interval <- fragment_size_interval
  metadata(gr)$bin_size <- bin_size
  metadata(gr)$window_size <- window_size 
  metadata(gr)$expand <- expand
  metadata(gr)$n_intervals <- n_intervals
  metadata(gr)$num_samples <- num_samples
  metadata(gr)$n_bins_per_window <- n_bins_per_window 
  metadata(gr)$step_size <- step_size

	gr

} # readFragmentSizeMatrix


makeData <- function(x, min_reads_per_window = 10, min_reads_coverage = 20){

	n_bins_per_window <- metadata(x)$window_size / metadata(x)$bin_size

	# the index of the 10-bp bins in [-200, +200] window
	starts <- seq(1, metadata(x)$expand / metadata(x)$bin_size - n_bins_per_window + 1, by = metadata(x)$step_size / metadata(x)$bin_size)
	ends <- starts + n_bins_per_window - 1
	n_block <- length(starts)

	# reads coverage of moving 400-bp window for each bin in [-200, +200] window (aka, core bins)
	# with step size 200 bp
	gr <- Reduce('c', lapply(1:n_block, function(j){

		# the index of the bins covering the [-200, +200] window
		m <- starts[j]:ends[j]
		xj <- x
		ranges(xj) <- shift(resize(ranges(xj), fix = 'start', width = metadata(x)$window_size), metadata(x)$step_size * (j - 1))
		mcols(xj)$coverage <- mcols(xj)$coverage[, m, drop = FALSE]
		mcols(xj)$min_coverage <- rowMins(mcols(xj)$coverage)
		mcols(xj)$max_coverage <- rowMaxs(mcols(xj)$coverage)
		mcols(xj)$mean_coverage <- rowMeans(mcols(xj)$coverage)

		xj <- xj[mcols(xj)$mean_coverage >= min_reads_coverage]

		mcols(xj)$coverage <- (mcols(xj)$coverage - mcols(xj)$min_coverage) / (mcols(xj)$max_coverage - mcols(xj)$min_coverage)
		mcols(xj)$coverage[is.na(mcols(xj)$coverage)] <- 0

		m2 <- rep(m, metadata(x)$n_interval) + (rep(1:metadata(x)$n_interval, each = n_bins_per_window) - 1) * metadata(x)$expand / metadata(x)$bin_size 
		mcols(xj)$counts <- mcols(xj)$counts[, m2, drop = FALSE]
		mcols(xj)$num_reads <- rowSums(mcols(xj)$counts)
		xj <- xj[mcols(xj)$num_reads >= min_reads_per_window]
		mcols(xj)$block <- j
		xj
	}))

	metadata(gr)$expand <- NULL
	metadata(gr)$step_size <- NULL
	metadata(gr)$min_reads_per_window <- min_reads_per_window
	metadata(gr)$n_bins_per_window <- n_bins_per_window
	gr

} # makeTrainData
