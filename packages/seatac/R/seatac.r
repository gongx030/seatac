#' seatac
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians
#' @import tensorflow
#' @import keras
#' @importFrom futile.logger flog.info
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom abind abind
NULL

#' seatac 
#'
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(filenames, which = NULL, genome, latent_dim = 10, window_size = 2000, bin_size = 20, fragment_size_range = c(70, 500), fragment_size_interval = 10, epochs = 5, batch_size = 256, sigma = 1){

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
  seqlevels(which, pruning.mode = 'coarse') <- seqlevels(gr)

  fs <- getFragmentSizeMatrix(filenames, which, window_size, bin_size, fragment_size_range)

  input_dim <- dim(fs$X)[2]
  feature_dim <- dim(fs$X)[3]
  model <- cae(input_dim = input_dim, feature_dim = feature_dim, latent_dim = latent_dim)
  model %>% fit(fs$X, fs$X, shuffle = TRUE, epochs = epochs, batch_size = batch_size, validation_split = 0.1)

  Xp <- model %>% predict(fs$X[1:20, , ])
  browser()

  i <- 4
  plot(colMeans(Xp[i, , ]), main = 'predicted'); plot(colMeans(fs$X[i, , ]), main = 'observed')
  plot(colMeans(Xp[1, , ]), main = 'predicted'); points(colMeans(Xp[2, , ]), col = 'red'); points(colMeans(Xp[3, , ]), col = 'blue')

  image(Xp[i, , ]); image(fs$X[i, , ])

  model %>% train_on_batch(fs$X, fs$X)
}
