#' predict.vae
#' 
predict.vae <- function(model, gr, type = 'nucleosome', chunk_size = 2^17, batch_size = 256){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	window_dim2 <- window_dim / num_samples
	S <- matrix(1:window_dim, window_dim2, num_samples)
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	starts <- seq(1, window_dim2, by = chunk_size)
	ends <- starts + chunk_size - 1
	ends[ends > window_dim2] <- window_dim2
	n_chunk <- length(starts)

	# latent representation of each unique window
	Z <- matrix(NA, window_dim2, model$latent_dim)

	for (b in 1:n_chunk){

		flog.info(sprintf('%d/%d', b, n_chunk))
		i <- starts[b]:ends[b]

		x_list <- lapply(1:num_samples, function(j){
			mcols(gr[S[i, j]])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L))
		})

		Z[i, ]  <- model$latent %>% predict(x_list, batch_size = batch_size, verbose = 1)
	}

	mcols(gr)$latent <- Z[rep(1:window_dim2, num_samples), ]

	gr

} # predict.vae


#' encoding
#' 
encoding <- function(model, gr, batch_size = 256){

	window_dim <- length(gr)
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	x <- mcols(gr)$counts %>%
		as.matrix() %>%
		array_reshape(c(window_dim, model$feature_dim, model$input_dim, 1L))
	
	y <- mcols(gr)$coverage %>%
		array_reshape(c(window_dim, model$window_size, 1L))

	mcols(gr)$latent  <- model$latent %>% predict(list(x, y), batch_size = batch_size, verbose = 1)

	gr

} # encoding


