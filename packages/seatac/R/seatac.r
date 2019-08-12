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
	epochs = 50, 
	batch_size = 256,
	validation_split = 0.1
){

	flog.info(sprintf('window size: %d', metadata(x)$window_size))

  window_dim <- length(x)
  feature_dim <- metadata(x)$n_intervals
  input_dim <- metadata(x)$n_bins_per_window

	flog.info(sprintf('total number of training windows(window_dim): %d', window_dim))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

	vplot_inputs <- layer_input(shape = c(input_dim, feature_dim, 1L))

	vplot_output <- vplot_inputs %>% 
	  layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = shape(2L, 2L),
			activation = 'relu'
		) %>%
	  layer_batch_normalization() %>%
		layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = shape(2L, 2L),
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = shape(2L, 2L),
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_flatten() %>%
		layer_dense(units = 16L, activation = 'relu') %>%
		layer_dropout(rate = 0.3)

	coverage_inputs <- layer_input(shape = c(input_dim, 1L))

	coverage_output <- coverage_inputs %>%
		layer_conv_1d(
			filters = 8L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_flatten() %>%
		layer_dense(units = 16L, activation = 'relu') %>%
		layer_dropout(rate = 0.3)

	outputs <- layer_concatenate(list(vplot_output, coverage_output)) %>%
		layer_dense(units = 16L, activation = 'relu') %>%
		layer_dropout(rate = 0.3) %>%
		layer_dense(units = 1L, activation = 'sigmoid')

	model <- keras_model(
		inputs = list(vplot_inputs, coverage_inputs),
		outputs = outputs
	)

	model %>% compile(
		optimizer = 'rmsprop',
		loss = 'binary_crossentropy',
		metrics = c('accuracy')
	)
	print(summary(model))

	vplot <- mcols(x)$counts %>%
		as.matrix() %>%
		array_reshape(c(window_dim, input_dim, feature_dim, 1))

	cvg <- mcols(x)$coverage %>%
		array_reshape(c(window_dim, input_dim, 1))

	label <- mcols(x)$label %>% as.numeric()

	model %>% fit(
		x = list(vplot, cvg),
		y = label,
		batch_size = batch_size,
		epochs = epochs,
		validation_split = validation_split
	)

	model

} # seatac


