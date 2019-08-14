#' seatac
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians rowMins rowMaxs
#' @import futile.logger 
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
#' @importFrom gplots colorpanel
#' @importFrom abind abind
#' @import tensorflow
#' @import keras
#' @importFrom reticulate array_reshape
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
	prior = 'gmm',
	batch_effect = FALSE,
	epochs = 50, 
	batch_size = 256, 
	beta = 1,
	min_reads_per_window = 10,
	min_reads_coverage = 20
){

	flog.info(sprintf('window size: %d', metadata(x)$window_size))
	flog.info(sprintf('step size: %d', metadata(x)$step_size))

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('# mixture components(n_components):%d', n_components))
	flog.info(sprintf('modeling batch effect(batch_effect): %s', batch_effect))
	flog.info(sprintf('latent prior model: %s', prior))
	flog.info(sprintf('minimum PE read pairs in each window for training(min_reads_per_window): %d', min_reads_per_window))
	flog.info(sprintf('minimum average read coverage per window for training(min_reads_coverage): %d', min_reads_coverage))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	x <- makeData(x, min_reads_per_window = min_reads_per_window, min_reads_coverage = min_reads_coverage)

	window_dim <- length(x)
	feature_dim <- metadata(x)$n_intervals
	input_dim <- metadata(x)$n_bins_per_window

	flog.info(sprintf('total number of training windows(window_dim): %d', window_dim))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

	model <- vae(
		input_dim = input_dim, 
		feature_dim = feature_dim, 
		latent_dim = latent_dim, 
		n_components = n_components, 
		num_samples = metadata(x)$num_samples,
		prior = prior
	)
	model %>% fit(x, epochs = epochs, beta = beta)
	model

} # seatac

