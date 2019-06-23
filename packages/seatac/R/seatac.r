#' seatac
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians
#' @importFrom futile.logger flog.info
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom abind abind
#' @import ggbio
NULL

#' seatac 
#'
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(filenames, which = NULL, genome, window_size = 2000, bin_size = 20, fragment_size_range = c(70, 500), fragment_size_interval = 10, epochs = 5, batch_size = 256, gpu = TRUE){

  if (missing(filenames))
    stop('filenames are missing')

  num_samples <- length(filenames)

  existed <- file.exists(filenames)
  if (any(!existed))
    stop(sprintf('%s do not exist', paste(filenames[!existed], collapse = ', ')))

  is_pe <- sapply(filenames, testPairedEndBam)
  if(any(!is_pe)){
    stop(paste(filenames[!s_pe], collapse = ', '),"are not paired-end files.")
  }

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

  fs <- getFragmentSizeMatrix(filenames, which, window_size, bin_size, fragment_size_range, fragment_size_interval)

  sample_dim <- dim(fs$X)[1]
  input_dim <- dim(fs$X)[2]
  feature_dim <- dim(fs$X)[3]

  model <- build_model(input_dim = input_dim, feature_dim = feature_dim, gpu = gpu)

  outputs <- list(
    abind(array(0, dim = c(sample_dim, 1, feature_dim)), fs$X[, 1:(input_dim - 1), ], along = 2),
    fs$X,
    abind(fs$X[, 2:input_dim, ], array(0, dim = c(sample_dim, 1, feature_dim)), along = 2)
   )
  model %>% fit(fs$X, outputs, epochs = epochs, batch_size = batch_size, validation_split = 0.1)

  Xp <- (model %>% predict(fs$X, batch_size = batch_size, verbose = 1))[[2]]
  Xp <- aperm(Xp, c(2, 1, 3)) 
  dim(Xp) <- c(prod(dim(Xp)[1:2]), dim(Xp)[3])

  X <- fs$X
  X <- aperm(X, c(2, 1, 3)) 
  dim(X) <- c(prod(dim(X)[1:2]), dim(X)[3])

  mcols(fs$bins)$counts <- X
  mcols(fs$bins)$predicted_counts <- Xp

  fs$bins

} # seatac


