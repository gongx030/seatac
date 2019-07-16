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
	expand = 2000
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

	flog.info(sprintf('resizing each peak to [-%d, +%d] region', window_size, window_size))
	windows <- resize(windows, fix = 'center', width = window_size)	# resize to window_size

	Y <- lapply(cvg, function(c) c[windows] %>% as.matrix())
	Y <- unlist(Y)
	dim(Y) <- c(length(windows), window_size, num_samples)

	bins <- unlist(slidingWindows(windows, width = bin_size, step = bin_size))

  # find the # PE reads per window per sample
  num_reads <- do.call('cbind', lapply(1:num_samples, function(i) countOverlaps(windows, x[mcols(x)$group == i])))

  x <- subsetByOverlaps(x, windows)
	wb <- cbind(rep(1:length(windows), n_bins_per_window), 1:(length(windows)* n_bins_per_window))	# windows ~ bins, sorted by window

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
  mcols(windows)$num_reads <- num_reads
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
