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

  for (i in 1:num_samples){

    if (!file.exists(filenames[i]))
      stop(sprintf('%s do not exist', filenames[i]))

    if(!testPairedEndBam(filenames[i]))
      stop(sprintf('%s is not paired-end file.', filenames[i]))

    x <- readGAlignmentPairs(filenames[i], param = param)

    if (length(x) > 0){

      x <- as(x, 'GRanges') # convert GAlignmentPairs into GRanges where two paired end reads are merged
      x$fragment_size <- width(x)
      x$fragment_size <- as.numeric(cut(x$fragment_size, breaks))
      x <- x[!is.na(x$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"
      x <- resize(x, width = 1, fix = 'center')

      S <- sparseMatrix(1:length(x), x$fragment_size, dims = c(length(x), length(breaks)))  # read center ~ fragment size
      mm <- as.matrix(findOverlaps(bins, x))	# bins ~ read center
      A <- as(sparseMatrix(mm[, 1], mm[, 2], dims = c(length(bins), length(x))), 'dgCMatrix') # bins ~ read center
      X <- A %*% S  # bins ~ fragment size
      X[X > 0] <- 1

			mm <- as.matrix(findOverlaps(windows, bins))	# windows ~ bins
      B <- as(sparseMatrix(mm[, 1], mm[, 2], dims = c(length(windows), length(bins))), 'dgCMatrix') # windows ~ bins
			include_window <- rowSums(B %*% X) >= min_reads_per_window # windows that have enough reads

			if (any(include_window)){

				mm <- summary(t(B[include_window, , drop = FALSE]))[, 1:2]	# bins, windows, order by windows
				X <- X[mm[, 1], ]	# only include the bins overlap with included windows (with enough reads)
				X <- as.matrix(X) # convert dgCMatrix to matrix

				# convert X into an array with bathc_size ~ n_bins_per_window ~ n_intervals
				dim(X) <- c(n_bins_per_window, sum(include_window), n_intervals)
				X <- aperm(X, c(1, 3, 2)) # n_bins_per_window ~ n_intervals ~ bathc_size 
				dim(X) <- c(n_bins_per_window * n_intervals, sum(include_window))
				X <- aperm(X, c(2, 1))
				X <- as(X, 'dgCMatrix')

				win_i <- windows[include_window]
				mcols(win_i)$group <- i
				mcols(win_i)$num_reads <- rowSums(X)
				mcols(win_i)$counts <- X

				if (is.null(win)){
					win <- win_i
				}else{
					win <- c(win, win_i)
				}
			}
    }
  }

	if (!is.null(win)){
	  metadata(win)$fragment_size_range  <- fragment_size_range
		metadata(win)$fragment_size_interval <- fragment_size_interval
		metadata(win)$bin_size <- bin_size
		metadata(win)$window_size <- window_size 
		metadata(win)$n_bins_per_window <- n_bins_per_window
		metadata(win)$n_intervals <- n_intervals
		win
	}

} # getFragmentSizeMatrix

