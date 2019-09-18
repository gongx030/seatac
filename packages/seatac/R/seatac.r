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
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(
	x,			# GenomicRanges object
	latent_dim = 2, 
	epochs = 50, 
	batch_size = 256,
	min_reads_per_window = 10,
	sequence_dim = 16L,
	type = 'vae'
){

	flog.info(sprintf('window size: %d', metadata(x)$window_size))
	flog.info(sprintf('step size: %d', metadata(x)$step_size))

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	feature_dim <- metadata(x)$n_intervals
	input_dim <- metadata(x)$n_bins_per_window
	num_samples <- metadata(x)$num_samples
	window_size <- metadata(x)$window_size

	flog.info(sprintf('total number of input windows: %d', length(x)))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

	flog.info(sprintf('min PE reads per window(min_reads_per_window): %d', min_reads_per_window))
	x <- x[mcols(x)$num_reads >= min_reads_per_window]
	flog.info(sprintf('total number of training windows: %d', length(x)))

	flog.info(sprintf('model(type): %s', type))

	if (type == 'vae'){
		model <- vae(
			input_dim = input_dim, 
			feature_dim = feature_dim, 
			latent_dim = latent_dim, 
			num_samples = num_samples
		)
	}else if (type == 'cvae'){
		flog.info(sprintf('sequence dimension(sequence_dim):%d', sequence_dim))
		model <- cvae(
			input_dim = input_dim, 
			feature_dim = feature_dim, 
			latent_dim = latent_dim, 
			num_samples = num_samples,
			window_size = window_size ,
			sequence_dim = sequence_dim
		)
	}
	model %>% fit(x, epochs = epochs)
	model

} # seatac

