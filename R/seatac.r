#' seatac
#'
#' @import tensorflow
#' @import keras
#' @import tfprobability
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians rowMins rowMaxs
#' @import futile.logger 
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam indexBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
#' @importFrom gplots colorpanel
#' @importFrom abind abind
#' @importFrom reticulate array_reshape
NULL

#' seatac 
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(
	x,			# GenomicRanges object
	latent_dim = 2L, 
	epochs = 50,
	batch_size = 128,
	steps_per_epoch = 50,
	type = 'vae_baseline',
	...
){

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))

	feature_dim <- metadata(x)$n_intervals
	num_samples <- metadata(x)$num_samples

	flog.info(sprintf('total number of input windows: %d', length(x)))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

	flog.info(sprintf('input window size:%d', metadata(x)$window_size))

	input_dim <- metadata(x)$n_bins_per_window
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))

	model <- build_model(
		type = type,
		input_dim = input_dim, 
		feature_dim = feature_dim, 
		latent_dim = latent_dim, 
		num_samples = num_samples,
		window_size = metadata(x)$window_size ,
		...
	)

	model %>% fit(x, batch_size = batch_size, epochs = epochs, steps_per_epoch = steps_per_epoch)

} # seatac

