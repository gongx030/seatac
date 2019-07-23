#' readFragmentSizeMatrix
#' @export
#'
readFragmentSizeMatrix <- function(
	x, 
	windows, 
	expand = 500,
	bin_size = 10,
	fragment_size_range = c(50, 670), 
	fragment_size_interval = 20
){

  num_samples <- metadata(x)$num_samples
  n_bins_per_window <- expand / bin_size
	breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)

	n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval + 1

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

  metadata(gr)$fragment_size_range  <- fragment_size_range
  metadata(gr)$fragment_size_interval <- fragment_size_interval
  metadata(gr)$bin_size <- bin_size
  metadata(gr)$expand <- expand
  metadata(gr)$n_intervals <- n_intervals
  metadata(gr)$num_samples <- num_samples

	gr

} # readFragmentSizeMatrix


makeData <- function(x, window_size = 320, min_reads_per_window = 5){

	expand <- metadata(x)$expand
	bin_size <- metadata(x)$bin_size
	n_bins_per_window <- window_size / bin_size

	# the index of the 10-bp bins in [-160, +160] window
	js <- (n_bins_per_window / 2):(expand / bin_size - n_bins_per_window / 2)

	# reads coverage of moving 320-bp window for each bin in [-160, +160] window (aka, core bins)
	gr <- Reduce('c', lapply(js, function(j){
		# the index of the bins covering the [-160, +160] window
		m <- (j - n_bins_per_window / 2 + 1):(j + n_bins_per_window / 2)
		xj <- x
		ranges(xj) <- IRanges(start = start(xj) + bin_size * (j - 1), width = window_size)
		mcols(xj)$coverage <- mcols(xj)$coverage[, m]
		mcols(xj)$min_coverage <- rowMin(mcols(xj)$coverage)
		mcols(xj)$max_coverage <- rowMax(mcols(xj)$coverage)
		mcols(xj)$coverage <- (mcols(xj)$coverage - mcols(xj)$min_coverage) / (mcols(xj)$max_coverage - mcols(xj)$min_coverage)
		mcols(xj)$is_peak <- mcols(xj)$coverage[, n_bins_per_window / 2] == 1
		mcols(xj)$is_valley <- mcols(xj)$coverage[, n_bins_per_window / 2] == 0

		if (train){
			is_peak_or_valley <- mcols(xj)$is_peak | mcols(xj)$is_valley
			xj <- xj[!is.na(is_peak_or_valley) & is_peak_or_valley]
		}

		X <- mcols(xj)$counts %>% as.matrix()
		dim(X) <- c(length(xj), expand / bin_size, metadata(x)$n_intervals)
		X <- X[, m, ]
		dim(X) <- c(length(xj), n_bins_per_window * metadata(x)$n_intervals)
		X <- as(X, 'dgCMatrix')
		mcols(xj)$counts <- X
		mcols(xj)$num_reads <- rowSums(X)
		xj <- xj[mcols(xj)$num_reads > min_reads_per_window]
		xj
	}))

	metadata(gr)$window_size <- window_size
	metadata(gr)$min_reads_per_window <- min_reads_per_window
	metadata(gr)$n_bins_per_window <- n_bins_per_window
	gr

} # makeTrainData
