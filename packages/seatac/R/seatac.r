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
	type = 'vae',
	latent_dim = 2L, 
	hidden_dim = 8L,
	epochs = 50,
	batch_size = 256,
	steps_per_epoch = 10,
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

	H <- metadata(gr)$n_intervals * metadata(gr)$n_bins_per_window
	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_bins_per_window, each = metadata(gr)$n_intervals), dims = c(H, metadata(gr)$n_bins_per_window))
	fragment_size <- mcols(gr)$counts %*% A %>% as.matrix()
	mcols(x)$fragment_size <- (fragment_size - rowMins(fragment_size)) / (rowMaxs(fragment_size) - rowMins(fragment_size))

	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_intervals, metadata(gr)$n_bins_per_window), dims = c(H, metadata(gr)$n_bins_per_window))
	position <- mcols(gr)$counts %*% A %>% as.matrix()
	mcols(x)$position <- (position - rowMins(position)) / (rowMaxs(position) - rowMins(position))

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

