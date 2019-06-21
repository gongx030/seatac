#' readFragmentSize
#' @export

readFragmentSize <- function(filename, bins, fragment_size_range = c(0, 750), fragment_size_interval = 10){

  if (missing(filename))
    stop('filename are missing')

  if (!file.exists(filename))
    stop(sprintf('%s do not exist', filename))

  if(!testPairedEndBam(filename))
    stop(sprintf('%s is not paired-end file.', filename))

  breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
  M <- length(breaks)

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isProperPair = TRUE
  )

  param <- ScanBamParam(which = reduce(bins), flag = flag)
  x <- readGAlignmentPairs(filename, param = param)
  if (length(x) > 0){

    x <- as(x, 'GRanges') # convert GAlignmentPairs into GRanges where two paired end reads are merged
    x$fragment_size <- width(x)
    x$fragment_size <- as.numeric(cut(x$fragment_size, breaks))
    x <- x[!is.na(x$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"
    x <- resize(x, width = 1, fix = 'center')
    B <- sparseMatrix(1:length(x), x$fragment_size, dims = c(length(x), length(breaks)))  # read center ~ fragment size
    mm <- as.matrix(findOverlaps(bins, x))
    A <- as(sparseMatrix(mm[, 1], mm[, 2], dims = c(length(bins), length(x))), 'dgCMatrix') # bins ~ read center
    X <- A %*% B  # bins ~ fragment size
    X[X > 0] <- 1
    X <- as.matrix(X)
  }else
    X <- matrix(0, dims = c(length(bins), length(x)))

  X
} # readFragmentSize


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

  X <- NULL
  pos <- NULL
  group <- NULL
  for (i in 1:num_samples){
    Xi <- readFragmentSize(filenames[i], bins, fragment_size_range, fragment_size_interval)

    # convert Xi into an array with bathc_size ~ n_bins_per_window ~ n_intervals
    dim(Xi) <- c(n_bins_per_window, n_windows, n_intervals)
    Xi <- aperm(Xi, c(2, 1, 3))
    include <- apply(Xi, 1, sum) > min_reads_per_window
    Xi <- Xi[include, , ]

    group <- c(group, rep(i, nrow(Xi)))
    pos <- c(pos, which(include))
    X <- abind(X, Xi, along = 1)
  }

  list(X = X, bins = windows[pos], group = group)

} # getFragmentSizeMatrix

