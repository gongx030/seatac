#' getFragmentSizeMatrix
#' @export

getFragmentSizeMatrix <- function(filenames, which, genome, window_size = 2000, bin_size = 20, step_size = 20, fragment_size_range = c(0, 500), fragment_size_interval = 10, min_reads_per_window = 20){

  num_samples <- length(filenames)
  n_bins_per_window <- window_size / bin_size
  n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval + 1

  if (missing(which) || length(which) > 1)
    stop('which must contain one genomic range')

  # genomic ranges covered by the BAM files
  gr <- Reduce('intersect', lapply(filenames, function(f){
    x <- idxstatsBam(f)
    GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
  }))
  seqlengths(seqinfo(gr)) <- width(gr)
  genome(seqinfo(gr)) <- providerVersion(genome)
  seqlevels(which, pruning.mode = 'coarse') <- seqlevels(gr)
  seqlevels(gr, pruning.mode = 'coarse') <- seqlevels(which)
  seqlengths(seqinfo(which)) <-  seqlengths(seqinfo(gr))
  genome(seqinfo(which)) <-  genome(seqinfo(gr))

  # every window has the same width
  windows <- slidingWindows(which, width = window_size, step = step_size)
  windows <- Reduce('c', windows)
  windows <- windows[width(windows) == window_size]
  seqlevels(windows, pruning.mode = 'coarse') <- seqlevels(which)
  n_windows <- length(windows)

  bins <- slidingWindows(reduce(windows), width = bin_size, step = bin_size)
  bins <- unlist(bins)

  breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isProperPair = TRUE
  )
  param <- ScanBamParam(which = reduce(bins), flag = flag)

  group_window <- NULL
	win <- NULL	

  x <- lapply(1:num_samples, function(i){

    if (!file.exists(filenames[i]))
      stop(sprintf('%s do not exist', filenames[i]))

    if(!testPairedEndBam(filenames[i]))
      stop(sprintf('%s is not paired-end file.', filenames[i]))

    ga <- readGAlignmentPairs(filenames[i], param = param)
    flog.info(sprintf('reading %s', filenames[i]))

    if (length(ga) > 0)
      ga <- as(ga, 'GRanges') # convert GAlignmentPairs into GRanges where two paired end reads are merged
      mcols(ga)$group <- i
      ga
  })

  # !!! need to check whetehr there is any NULL
  x <- Reduce('c', x)
  x$fragment_size <- width(x)
  x$fragment_size <- as.numeric(cut(x$fragment_size, breaks))
  x <- x[!is.na(x$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"
  x <- resize(x, width = 1, fix = 'center') # the center genomic coordinate between two PE reads

  # find the # PE reads per window per sample
  num_reads <- do.call('cbind', lapply(1:num_samples, function(i) countOverlaps(windows, x[mcols(x)$group == i])))

  # each window must have at least min_reads_per_window reads in all samples
  include_window <- rowSums(num_reads >= min_reads_per_window) == num_samples

  if (any(include_window))
    flog.info(sprintf('selecting %d windows that have at least %d PE reads in all samples', sum(include_window), min_reads_per_window))
  else
    stop(sprintf('there is no window that have %d PE reads in all samples. Try smaller min_reads_per_window', sum(include_window), min_reads_per_window))

  num_reads <- num_reads[include_window, , drop = FALSE]
  windows <- windows[include_window]
  x <- subsetByOverlaps(x, windows)
  bins <- subsetByOverlaps(bins, windows)
  wb <- as.matrix(findOverlaps(windows, bins))	# windows ~ bins, sorted by window

  X <- Reduce('cbind', lapply(1:num_samples, function(i){
    xi <- x[mcols(x)$group == i]  # reads from group i
    CF <- sparseMatrix(i = 1:length(xi), j = xi$fragment_size, dims = c(length(xi), length(breaks)))  # read center ~ fragment size
    BC <- as.matrix(findOverlaps(bins, xi))	# bins ~ read center
    BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(xi))), 'dgCMatrix') # bins ~ read center
    BF <- BC %*% CF  # bins ~ fragment size
    BF[BF > 0] <- 1
    BF <- as.matrix(BF[wb[, 2], ])

    # convert BF into an array with bathc_size ~ n_bins_per_window ~ n_intervals
    dim(BF) <- c(n_bins_per_window, sum(include_window), n_intervals)
    BF <- aperm(BF, c(1, 3, 2)) # n_bins_per_window ~ n_intervals ~ bathc_size 
    dim(BF) <- c(n_bins_per_window * n_intervals, sum(include_window))
    BF <- aperm(BF, c(2, 1))  # bathc_size ~ n_bins_per_window * n_intervals
    BF <- as(BF, 'dgCMatrix')
    BF
  }))

  mcols(windows)$counts <- X
  mcols(windows)$num_reads <- num_reads
  metadata(windows)$fragment_size_range  <- fragment_size_range
  metadata(windows)$fragment_size_interval <- fragment_size_interval
  metadata(windows)$bin_size <- bin_size
  metadata(windows)$window_size <- window_size 
  metadata(windows)$n_bins_per_window <- n_bins_per_window
  metadata(windows)$n_intervals <- n_intervals
  windows

} # getFragmentSizeMatrix


readFragmentSize <- function(
	filenames, 
	which = NULL, 
	genome, 
	window_size = 320, 
	bin_size = 10, 
	fragment_size_range = c(50, 670), 
	fragment_size_interval = 20, 
	min_reads_per_window = 50
){

	validate_bam(filenames)

  num_samples <- length(filenames)

  input_dim <- window_size / bin_size
	flog.info(sprintf('window size(window_size): %d', window_size))
	flog.info(sprintf('bin size(bin_size): %d', bin_size))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))

  gr <- getFragmentSizeMatrix(filenames, which, genome, window_size, bin_size, bin_size, fragment_size_range, fragment_size_interval, min_reads_per_window)
  metadata(gr)$num_samples <- num_samples
  gr
} # readFragmentSize
