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
	epochs = 50, 
	batch_size = 256,
	min_reads_per_window = 5,
	min_reads_coverage = 5
){

	flog.info(sprintf('window size: %d', metadata(x)$window_size))
	flog.info(sprintf('step size: %d', metadata(x)$step_size))

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	window_dim <- length(x)
	feature_dim <- metadata(x)$n_intervals
	input_dim <- metadata(x)$n_bins_per_window
	window_size <- metadata(x)$window_size

	flog.info(sprintf('total number of training windows(window_dim): %d', window_dim))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))
	flog.info(sprintf('min PE reads per window(min_reads_per_window): %d', min_reads_per_window))
	flog.info(sprintf('min reads coverage per window(min_reads_coverage): %d', min_reads_coverage))

	model <- vae(
		input_dim = input_dim, 
		feature_dim = feature_dim, 
		latent_dim = latent_dim, 
		window_size = window_size,
		num_samples = metadata(x)$num_samples
	)

	train <- mcols(windows)$num_reads >= min_reads_per_window & mcols(windows)$mean_coverage >= min_reads_coverage
	model %>% fit(x[train], epochs = epochs)
	model

} # seatac

