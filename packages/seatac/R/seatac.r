#' seatac
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians
#' @importFrom futile.logger flog.info
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
#' @importFrom gplots colorpanel
#' @importFrom abind abind
#' @import tensorflow
#' @import keras
#' @import tfprobability
NULL


#' seatac 
#'
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(
	x,			# GenomicRanges object
	latent_dim = 2, 
	n_components = 10, 
	window_size = 320,
	prior = 'gmm',
	batch_effect = FALSE,
	epochs = 50, 
	batch_size = 256, 
	beta = 1
){

	flog.info('preparing training data')
	train <- makeData(x, window_size = window_size, train = TRUE)

  window_dim <- length(train)
  feature_dim <- metadata(train)$n_intervals
  input_dim <- metadata(train)$n_bins_per_window

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('# mixture components(n_components):%d', n_components))
	flog.info(sprintf('total number of training windows(window_dim): %d', window_dim))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))
	flog.info(sprintf('modeling batch effect(batch_effect): %s', batch_effect))
	flog.info(sprintf('latent prior model: %s', prior))

  model <- vae(
    input_dim = input_dim, 
    feature_dim = feature_dim, 
    latent_dim = latent_dim, 
    n_components = n_components, 
    num_samples = metadata(x)$num_samples,
		prior = prior
  )

	model %>% fit(train, epochs = epochs, beta = beta)

	mcols(x)$latent <- model %>% predict(x, window_size = window_size)

	metadata(x)$model <- model
  x	

} # seatac


