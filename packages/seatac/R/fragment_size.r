#' getFragmentSizeMatrix
#' @export

getFragmentSizeMatrix <- function(filenames, which, window_size = 2000, bin_size = 20, fragment_size_range = c(0, 500), fragment_size_interval = 10, min_reads_per_window = 20){

  num_samples <- length(filenames)
  n_bins_per_window <- window_size / bin_size
  n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval + 1

  if (missing(which) || length(which) > 1)
    stop('which must contain one genomic range')

  # every window has the same width
  windows <- slidingWindows(which, width = window_size, step = window_size)
  windows <- Reduce('c', windows)
  windows <- windows[width(windows) == window_size]
  seqlevels(windows, pruning.mode = 'coarse') <- seqlevels(which)
  n_windows <- length(windows)

  bins <- slidingWindows(windows, width = bin_size, step = bin_size)
  bins <- unlist(bins)

  breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isProperPair = TRUE
  )
  param <- ScanBamParam(which = reduce(bins), flag = flag)

  X <- NULL
  group_window <- NULL
  group_bin <- NULL
  cvg <- NULL

  for (i in 1:num_samples){

    if (!file.exists(filenames[i]))
      stop(sprintf('%s do not exist', filenames[i]))

    if(!testPairedEndBam(filenames[i]))
      stop(sprintf('%s is not paired-end file.', filenames[i]))

    x <- readGAlignmentPairs(filenames[i], param = param)
    cvg_i <- sum(coverage(x)[bins])  # total reads coverage at each bin
    names(cvg_i) <- NULL
    cvg <- c(cvg, cvg_i)

    if (length(x) > 0){
      x <- as(x, 'GRanges') # convert GAlignmentPairs into GRanges where two paired end reads are merged
      x$fragment_size <- width(x)
      x$fragment_size <- as.numeric(cut(x$fragment_size, breaks))
      x <- x[!is.na(x$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"
      x <- resize(x, width = 1, fix = 'center')
      S <- sparseMatrix(1:length(x), x$fragment_size, dims = c(length(x), length(breaks)))  # read center ~ fragment size
      mm <- as.matrix(findOverlaps(bins, x))
      A <- as(sparseMatrix(mm[, 1], mm[, 2], dims = c(length(bins), length(x))), 'dgCMatrix') # bins ~ read center
      Xi <- A %*% S  # bins ~ fragment size
      Xi[Xi > 0] <- 1
      Xi <- as.matrix(Xi)
    }else{
      Xi <- matrix(0, dims = c(length(bins), length(x)))
    }

    # convert Xi into an array with bathc_size ~ n_bins_per_window ~ n_intervals
    dim(Xi) <- c(n_bins_per_window, n_windows, n_intervals)
    Xi <- aperm(Xi, c(2, 1, 3))

    #include <- apply(Xi, 1, sum) > min_reads_per_window # selecting windows with enough reads
    #Xi <- Xi[include, , ]
    #Bi <- Bi[include, ]

    group_window <- c(group_window, rep(i, nrow(Xi)))
    group_bin <- c(group_bin, rep(i, n_bins_per_window * n_windows))
    X <- abind(X, Xi, along = 1)
  }

  bins <- rep(bins, num_samples)
  metadata(bins)$fragment_size_range  <- fragment_size_range
  metadata(bins)$fragment_size_interval <- fragment_size_interval
  metadata(bins)$bin_size <- bin_size
  mcols(bins)$coverage <- cvg   # coverage of each bin: bins ~ group
  mcols(bins)$groups <- group_bin   # group index for each bin

  list(
    X = X,  # input data: samples ~ input_dim ~ feature_dim
    group_window = group_window,  # group index for each window (used for training)
    bins = bins # GRange for each bin
  )

} # getFragmentSizeMatrix

