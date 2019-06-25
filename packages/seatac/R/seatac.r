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
#' @import ggbio
NULL

#' seatac 
#'
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(filenames, which = NULL, genome, latent_dim = 10, n_components = 5, window_size = 2000, bin_size = 20, fragment_size_range = c(70, 500), fragment_size_interval = 10, min_reads_per_window = 20, epochs = 5, batch_size = 256, steps_per_epoch = 10){

  if (missing(filenames))
    stop('filenames are missing')

  num_samples <- length(filenames)

  existed <- file.exists(filenames)
  if (any(!existed))
    stop(sprintf('%s do not exist', paste(filenames[!existed], collapse = ', ')))

  is_pe <- sapply(filenames, testPairedEndBam)
  if(any(!is_pe)){
    stop(paste(filenames[!s_pe], collapse = ', '),"are not paired-end files.")
  }

  # genomic ranges covered by the BAM files
  gr <- Reduce('intersect', lapply(filenames, function(f){
    x <- idxstatsBam(f)
    GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
  }))
  seqlengths(seqinfo(gr)) <- width(gr)
  genome(seqinfo(gr)) <- providerVersion(genome)
  seqlevels(which, pruning.mode = 'coarse') <- seqlevels(gr)
  seqlevels(gr, pruning.mode = 'coarse') <- seqlevels(which)
  seqlengths(seqinfo(which)) <-  seqlengths(seqinfo(gr))
  genome(seqinfo(which)) <-  genome(seqinfo(gr))

	flog.trace(sprintf('latent dimension(latent_dim):%d', latent_dim))
	flog.trace(sprintf('# mixture components(n_components):%d', n_components))

  fs <- getFragmentSizeMatrix(filenames, which, window_size, bin_size, fragment_size_range, fragment_size_interval, min_reads_per_window)

  sample_dim <- dim(fs$X)[1]
  input_dim <- dim(fs$X)[2]
  feature_dim <- dim(fs$X)[3]

	flog.trace(sprintf('total number of training windows(sample_dim): %d', sample_dim))
	flog.trace(sprintf('# bins per window(input_dim): %d', input_dim))
	flog.trace(sprintf('# features per bin(feature_dim): %d', feature_dim))

  model <- gmm_vae(input_dim = input_dim, feature_dim = feature_dim, latent_dim = latent_dim, n_components = n_components)

  train_dataset <- fs$X %>%
		k_expand_dims() %>%
		tf$cast(tf$float32) %>%
		tensor_slices_dataset() %>%
		dataset_shuffle(buffer_size = batch_size * steps_per_epoch) %>%
		dataset_batch(batch_size)

	model %>% fit_vae(train_dataset, epochs = epochs, steps_per_epoch = steps_per_epoch)

	Z <- model$encoder(fs$X %>% k_expand_dims() %>% tf$cast(tf$float32))$loc
	P <- model$latent_prior_model(NULL)$components_distribution$prob(Z %>% tf$reshape(shape(sample_dim, 1, latent_dim))) %>% as.matrix()
	Xp <- model$decoder(Z)$distribution$probs %>% tf$reshape(shape(sample_dim, input_dim, feature_dim)) %>% as.array()


	browser()
	par(mfrow = c(4, 6))
	lapply(1:n_components, function(k) image(colSums(Xp[max.col(P) == k, , ], dims = 1), main = k))
	lapply(1:n_components, function(k) image(colSums(fs$X[max.col(P) == k, , ], dims = 1), main = k))

	i <- 1
	par(mfrow = c(4, 6)); lapply(c(1, 50, 100, 200, 300, 500, 1000), function(i){image(fs$X[i, , ]); image(Xp[i, , ])})


	g <- array(0, dim = c(sample_dim, 1, feature_dim))
  outputs <- list(
    abind(g, fs$X[, 1:(input_dim - 1), ], along = 2),
    fs$X,
    abind(fs$X[, 2:input_dim, ], g, along = 2)
   )
  model %>% fit(fs$X, outputs, epochs = epochs, batch_size = batch_size, validation_split = 0.1, shuffle = TRUE)

} # seatac


