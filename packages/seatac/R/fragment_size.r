#' readFragmentSize
#' @export

readFragmentSize <- function(
	filenames, 
	peaks,
	genome, 
	window_size = 320,
	bin_size = 10, 
	fragment_size_range = c(50, 670), 
	fragment_size_interval = 20, 
	expand = 2000,
	min_reads_per_window = 5
){

	validate_bam(filenames)

  num_samples <- length(filenames)
  input_dim <- window_size / bin_size
	flog.info(sprintf('window size(window_size): %d', window_size))
	flog.info(sprintf('bin size(bin_size): %d', bin_size))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))

  n_bins_per_window <- window_size / bin_size
  n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval + 1

  # genomic ranges covered by the BAM files
  gr <- Reduce('intersect', lapply(filenames, function(f){
    x <- idxstatsBam(f)
    GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
  }))
  seqlengths(seqinfo(gr)) <- width(gr)
  genome(seqinfo(gr)) <- providerVersion(genome)

	# set the seqlevels and seqlengths from BAM files to "which"
  seqlevels(peaks, pruning.mode = 'coarse') <- seqlevels(gr)
  seqlevels(gr, pruning.mode = 'coarse') <- seqlevels(peaks)
  seqlengths(seqinfo(peaks)) <-  seqlengths(seqinfo(gr))
  genome(seqinfo(peaks)) <-  genome(seqinfo(gr))

  breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)

	flog.info(sprintf('expanding each peak to [-%d, +%d] region', expand / 2, expand / 2))
	windows <- resize(peaks, fix = 'center', width = expand)	# expand the input peaks

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isProperPair = TRUE
  )
  param <- ScanBamParam(which = reduce(windows), flag = flag, what = 'isize')

	# Read the PE reads overlapping with specified windows
  x <- lapply(1:num_samples, function(i){
    flog.info(sprintf('reading %s', filenames[i]))
    readGAlignments(filenames[i], param = param)
	})

	cvg <- lapply(x, coverage)	# reads coverage

	x <- lapply(1:num_samples, function(i){
		ga <- x[[i]]
		ga <- ga[strand(ga) == '+']
		# compute the center point between PE reads
		# this is faster than using GAlignmentPairs
		GRanges(
			seqnames = seqnames(ga), 
			range = IRanges(start(ga) + round(mcols(ga)$isize / 2), width = 1), 
			isize = mcols(ga)$isize, 
			group = i
		)
	})

  # !!! need to check whetehr there is any NULL
  x <- Reduce('c', x)
  x$fragment_size <- as.numeric(cut(x$isize, breaks))	# discretize the fragment size
  x <- x[!is.na(x$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"

	# find the peaks and valleys within each [-160, +160] window
	# For each window, we find a peak and a valley
	flog.info(sprintf('finding the peaks and valleys in each [-%d, +%d] region centered at the input summits', window_size / 2, window_size / 2))

	bins <- windows %>%
		resize(fix = 'center', width = window_size * 2 ) %>%
		slidingWindows(width = bin_size, step = bin_size) %>%
		unlist()

	mcols(bins)$peak_id <- rep(1:length(windows), each = n_bins_per_window * 2)
	mcols(bins)$bin_id <- rep(1:(n_bins_per_window * 2), length(windows))

	# coverage of each 10-bp bin in [-320, +320] window
	Y <- lapply(cvg, function(c) matrix(mean(c[bins]), nrow = length(windows), ncol = n_bins_per_window * 2, byrow = TRUE))

	# the index of the bins in [-160, +160] window
	js <- (n_bins_per_window / 2): (2 * n_bins_per_window - n_bins_per_window / 2)

	# bin index in Y
	B <- matrix(1:length(bins), length(windows), ncol = n_bins_per_window * 2, byrow = TRUE)

	# reads coverage of moving 320-bp window for each bin in [-160, +160] window (aka, core bins)
	Y <- lapply(1:num_samples, function(i) 
		do.call('rbind', lapply(js, function(j){
			m <- (j - n_bins_per_window / 2 + 1):(j + n_bins_per_window / 2)
			Y[[i]][, m]
		}))
	)
	max_coverage <- do.call('cbind', lapply(Y, rowMax))
	min_coverage <- do.call('cbind', lapply(Y, rowMin))

	bins <- bins[Reduce('c', lapply(js, function(j) B[, j]))]

	# find the windows that have reads coverage
	window_with_high_coverage <- rowSums(max_coverage > min_coverage) == num_samples
	max_coverage <- max_coverage[window_with_high_coverage, , drop = FALSE]
	min_coverage <- min_coverage[window_with_high_coverage, , drop = FALSE]
	Y <- lapply(1:num_samples, function(i) Y[[i]][window_with_high_coverage, , drop = FALSE])
	bins <- bins[window_with_high_coverage]

	# for each row, scale to [0, 1]
	Y <- lapply(1:num_samples, function(i) (Y[[i]] - min_coverage[, i]) / (max_coverage[, i] - min_coverage[, i]))

	# determine the peak/valley based on the scaled coverage density
	is_peak <- do.call('cbind', lapply(1:num_samples, function(i) Y[[i]][, n_bins_per_window / 2] == 0))
	is_valley <- do.call('cbind', lapply(1:num_samples, function(i) Y[[i]][, n_bins_per_window / 2] == 1))

	# include the windows that is a peak or a valley in at least one sample
	is_peak_or_valley <- rowSums(is_peak) > 0 | rowSums(is_valley) > 0

	# include only the peak/valley windows
	bins <- bins[is_peak_or_valley]
	is_peak <- is_peak[is_peak_or_valley, , drop = FALSE]
	is_valley <- is_valley[is_peak_or_valley, , drop = FALSE]
	max_coverage <- max_coverage[is_peak_or_valley, , drop = FALSE]
	min_coverage <- min_coverage[is_peak_or_valley, , drop = FALSE]
	Y <- lapply(1:num_samples, function(i) Y[[i]][is_peak_or_valley, , drop = FALSE])
	Y <- unlist(Y)
	dim(Y) <- c(length(bins), n_bins_per_window, num_samples)

	windows <- resize(bins, fix = 'center', width = window_size)	# resize to window_size
	bins <- unlist(slidingWindows(windows, width = bin_size, step = bin_size))  # all bins within each 320 bp window

  # find the # PE read centers per window per sample
  x <- subsetByOverlaps(x, windows)
	wb <- cbind(rep(1:length(windows), n_bins_per_window), 1:(length(windows)* n_bins_per_window))	# windows ~ bins, sorted by window
	num_reads <- do.call('cbind', lapply(1:num_samples, function(i) countOverlaps(windows, x[mcols(x)$group == i])))

  X <- Reduce('cbind', lapply(1:num_samples, function(i){
    xi <- x[mcols(x)$group == i]  # reads from group i
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

  mcols(windows)$counts <- X
  mcols(windows)$coverage <- Y
  mcols(windows)$max_coverage <- max_coverage 
  mcols(windows)$min_coverage <- min_coverage 
  mcols(windows)$is_peak <- is_peak
  mcols(windows)$is_valley <- is_valley 
  mcols(windows)$num_reads <- num_reads

	windows <- windows[rowSums(mcols(windows)$num_reads > min_reads_per_window) == num_samples]

  metadata(windows)$fragment_size_range  <- fragment_size_range
  metadata(windows)$fragment_size_interval <- fragment_size_interval
  metadata(windows)$bin_size <- bin_size
  metadata(windows)$window_size <- window_size 
  metadata(windows)$n_bins_per_window <- n_bins_per_window
  metadata(windows)$n_intervals <- n_intervals
  metadata(windows)$num_samples <- num_samples
	metadata(windows)$expand <- expand
  windows

} # readFragmentSize
