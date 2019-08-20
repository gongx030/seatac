#' predict.vae
#' 
predict.vae <- function(model, gr, type = 'nucleosome', batch_size = 256){

	window_dim <- length(gr)
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	x <- mcols(gr)$counts %>%
		as.matrix() %>%
		array_reshape(c(window_dim, model$input_dim, model$feature_dim, 1L))

	y <- mcols(gr)$coverage %>%
		array_reshape(c(window_dim, model$window_size, 1L))

	r <- model$nucleosome %>% predict(list(x, y), batch_size = batch_size, verbose = 1)

	mcols(gr)$fitted_counts <- r[[1]] %>% array_reshape(dim = c(window_dim, model$input_dim *  model$feature_dim), 'C')
	mcols(gr)$fitted_coverage <- r[[2]][, , 1]
	mcols(gr)$fitted_nucleosome_score <- r[[3]]

	gr

} # predict.vae


