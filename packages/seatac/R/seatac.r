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
	filenames, 
	which = NULL, 
	genome, 
	latent_dim = 10, 
	n_components = 5, 
	window_size = 320, 
	bin_size = 10, 
	fragment_size_range = c(50, 670), 
	fragment_size_interval = 20, 
	min_reads_per_window = 50, 
	epochs = 5, 
	batch_size = 256, 
	steps_per_epoch = 10
){

	validate_bam(filenames)

  num_samples <- length(filenames)

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('# mixture components(n_components):%d', n_components))

  win <- getFragmentSizeMatrix(filenames, which, genome, window_size, bin_size, bin_size, fragment_size_range, fragment_size_interval, min_reads_per_window)
  sample_dim <- length(win)
  input_dim <- window_size / bin_size
  feature_dim <- metadata(win)$n_bins_per_window

	flog.info(sprintf('total number of training windows(sample_dim): %d', sample_dim))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

  model <- gmm_vae(input_dim = input_dim, feature_dim = feature_dim, latent_dim = latent_dim, n_components = n_components, num_samples = num_samples)
	model %>% fit(win, epochs = epochs, steps_per_epoch = steps_per_epoch, batch_size = batch_size)

	x <- mcols(win)$counts %>%
		as.matrix() %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(sample_dim, model$input_dim, model$feature_dim)) %>%
		tf$expand_dims(axis = 3L)

	g <- mcols(win)$group - 1 %>%   # group index of current batch
		tf$cast(tf$int32)

	# latent representation
	Z <- model$encoder(x)$loc
	P <- model$latent_prior_model(NULL)$components_distribution$log_prob(
		Z %>% tf$reshape(shape(sample_dim, 1, model$latent_dim))
	) %>% as.matrix()
	mcols(win)$cluster <- max.col(P)
	metadata(win)$model <- model
	win

} # seatac


