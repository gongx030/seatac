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
#' @import Gviz
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
	min_reads_per_window_train = 50, 
	min_reads_per_window_predict = 0, 
	epochs = 5, 
	batch_size = 256, 
	steps_per_epoch = 10,
	num_windows_per_block = 1000 # size of genomic windows for prediction per chunk
){

	validate_bam(filenames)

  num_samples <- length(filenames)

	flog.info(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.info(sprintf('# mixture components(n_components):%d', n_components))

  win <- getFragmentSizeMatrix(filenames, which, genome, window_size, bin_size, bin_size, fragment_size_range, fragment_size_interval, min_reads_per_window_train)

  sample_dim <- dim(mcols(win)$counts)[1]
  input_dim <- dim(mcols(win)$counts)[2]
  feature_dim <- dim(mcols(win)$counts)[3]

	flog.info(sprintf('total number of training windows(sample_dim): %d', sample_dim))
	flog.info(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.info(sprintf('# features per bin(feature_dim): %d', feature_dim))

  model <- gmm_vae(input_dim = input_dim, feature_dim = feature_dim, latent_dim = latent_dim, n_components = n_components)

  train_dataset <- mcols(win)$counts %>%
		k_expand_dims() %>%
		tf$cast(tf$float32) %>%
		tensor_slices_dataset() %>%
		dataset_shuffle(buffer_size = batch_size * steps_per_epoch) %>%
		dataset_batch(batch_size)

	model %>% fit(train_dataset, epochs = epochs, steps_per_epoch = steps_per_epoch)

	model$which <- which
	model$genome <- genome
	model$window_size <- window_size
	model$bin_size <- bin_size
	model$fragment_size_range <- fragment_size_range
	model$fragment_size_interval <- fragment_size_interval
	model$min_reads_per_window <- min_reads_per_window_train 

	gr <- model %>% predict(filenames, which, num_windows_per_block = num_windows_per_block, min_reads_per_window = min_reads_per_window_predict)	
	gr

} # seatac


