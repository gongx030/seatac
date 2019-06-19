#' readFragmentSize
#' @export

readFragmentSize <- function(filename, bins, fragment_size_range = c(70, 750), fragment_size_interval = 10){

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

  }else
    X <- Matrix(0, dims = c(length(bins), length(x)))

  X
} # readFragmentSize


getFragmentSizeMatrix <- function(filenames, bins, expand = 2000, fragment_size_range = c(70, 500), fragment_size_interval = 10, min_reads_per_bin = 20, library_size = 1e5){

  num_samples <- length(filenames)

  bins2 <- resize(bins, expand, fix = 'center')
  bins2 <- trim(bins2)

  X <- NULL
  pos <- NULL
  group <- NULL
  for (i in 1:num_samples){
    Xi <- readFragmentSize(filenames[i], bins2, fragment_size_range, fragment_size_interval)
    included <- rowSums(Xi) >= min_reads_per_bin
    Xi <- Xi[included, ]

#    Xi <- Xi * library_size / sum(Xi)

    group <- c(group, rep(i, nrow(Xi)))
    pos <- c(pos, which(included))
    X <- rbind(X, Xi)
  }
  list(X = X, bins = bins[pos], group = group)

} # getFragmentSizeMatrix

