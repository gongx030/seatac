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
	epochs = 5, 
	batch_size = 256, 
	steps_per_epoch = 10
){

  sample_dim <- length(x)
  feature_dim <- metadata(x)$n_bins_per_window
  input_dim <- metadata(x)$window_size / metadata(x)$bin_size

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('# mixture components(n_components):%d', n_components))
	flog.info(sprintf('total number of training windows(sample_dim): %d', sample_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

	flog.info(sprintf('prior: %s', prior))

  if (prior == 'gmm'){
    model <- gmm_vae(
      input_dim = input_dim, 
      feature_dim = feature_dim, 
      latent_dim = latent_dim, 
      n_components = n_components, 
      num_samples = metadata(x)$num_samples
    )
  }

	model %>% fit(x, epochs = epochs, steps_per_epoch = steps_per_epoch, batch_size = batch_size)

	X <- mcols(x)$counts %>%
		as.matrix() %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(sample_dim, model$input_dim, model$feature_dim)) %>%
		tf$expand_dims(axis = 3L)

	g <- mcols(x)$group - 1 %>%   # group index of current batch
		tf$cast(tf$int32)

	# latent representation
	Z <- model$encoder(X)$loc
	P <- model$latent_prior_model(NULL)$components_distribution$log_prob(
		Z %>% tf$reshape(shape(sample_dim, 1, model$latent_dim))
	) %>% as.matrix()
	mcols(x)$cluster <- max.col(P)
	metadata(x)$model <- model
  x	

} # seatac


