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
#' @importFrom gplots colorpanel
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
  x,
	latent_dim = 10, 
	n_components = 5, 
  prior = 'gmm',
	batch_effect = FALSE,
	epochs = 5, 
	batch_size = 256, 
	steps_per_epoch = 10
){

  window_dim <- length(x)
  feature_dim <- metadata(x)$n_bins_per_window
  input_dim <- metadata(x)$window_size / metadata(x)$bin_size

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('# mixture components(n_components):%d', n_components))
	flog.info(sprintf('total number of training windows(window_dim): %d', window_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))
	flog.info(sprintf('modeling batch effect(batch_effect): %s', batch_effect))
	flog.info(sprintf('latent prior model: %s', prior))

  model <- vae(
    input_dim = input_dim, 
    feature_dim = feature_dim, 
    latent_dim = latent_dim, 
    n_components = n_components, 
    num_samples = metadata(x)$num_samples,
		batch_effect = batch_effect,
		prior = prior
  )
	browser()

	model %>% fit(x, epochs = epochs, steps_per_epoch = steps_per_epoch, batch_size = batch_size)
	mcols(x)$cluster <- model %>% predict(x, batch_size = batch_size)
	metadata(x)$model <- model
  x	

} # seatac


