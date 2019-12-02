#' predict.vae
#' 
predict.vae <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z
#	nucleosome <- matrix(NA, window_dim, model$input_dim) 	

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		z[i, ] <- model$encoder(x)$mean() %>% as.matrix()
#		nucleosome[i, ] <- model$decoder(z[i, ])[[1]] %>% as.matrix()

	}

	mcols(gr)$latent <- z
#	mcols(gr)$nucleosome <- nucleosome

	gr

} # predict.vae
