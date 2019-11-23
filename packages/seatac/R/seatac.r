#' seatac
#'
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
#' @import tensorflow
#' @import keras
#' @importFrom reticulate array_reshape
#' @importFrom DECIPHER AlignSeqs
#' @importFrom TFBSTools consensusMatrix seqLogo
NULL

#' seatac 
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(
	x,			# GenomicRanges object
	latent_dim = 2, 
	window_size = 640,
	min_reads_per_window = 20,
	epochs = 50,
	batch_size = 128,
	steps_per_epoch = 20,
	sequence_dim = 32,
	n_components = 15,
	core_width = 100,
	batch_size_prediction = 64,
	genome,
	...
){

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))

	feature_dim <- metadata(x)$n_intervals
	num_samples <- metadata(x)$num_samples

	flog.info(sprintf('total number of input windows: %d', length(x)))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

	flog.info(sprintf('input window size:%d', metadata(x)$window_size))
	flog.info(sprintf('training window size:%d', window_size))

	train <- sample_windows(x, window_size, min_reads_per_window, epochs = epochs, batch_size = batch_size, steps_per_epoch = steps_per_epoch)
	mcols(train)$sequence <- getSeq(genome, train)

	input_dim <- metadata(train)$n_bins_per_window
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))

	model <- gmm_cvae(
		input_dim = input_dim, 
		feature_dim = feature_dim, 
		latent_dim = latent_dim, 
		num_samples = num_samples,
		window_size = window_size ,
		sequence_dim = sequence_dim,
		n_components = n_components
	)

	model %>% fit(train)
	train <- model %>% predict(train, batch_size = batch_size_prediction)
	model$input_data <- x
	model$train_data <- train
	model

} # seatac

