#' predict.vae_baseline
#'
predict.vae_baseline <- function(model, gr, batch_size = 128){

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

	mcols(gr)$latent <- matrix(NA, window_dim, model$latent_dim)
	mcols(gr)$fragment_size_profile <- matrix(NA, window_dim, model$feature_dim)
	mcols(gr)$nfr_reads <- matrix(NA, window_dim, model$input_dim)
	mcols(gr)$mono_nucleosome_reads <- matrix(NA, window_dim, model$input_dim)
	is_nfr <- which(metadata(gr)$nfr)
	is_mono <- which(metadata(gr)$mono_nucleosome)

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		z <- x %>% model$encoder()
		z <- z$mean()

		y <- z %>% model$decoder()
		y <- y$mean()

		mcols(gr)$latent[i, ] <- z %>% as.matrix()	# latent space
		mcols(gr)$fragment_size_profile[i, ] <- y %>% tf$reduce_sum(axis = c(2L, 3L)) %>% as.matrix()	# fragment size profile
		mcols(gr)$nfr_reads[i, ] <- y[, is_nfr, , ] %>% tf$reduce_sum(axis = c(1L, 3L)) %>% as.matrix()
		mcols(gr)$mono_nucleosome_reads[i, ] <- y[, is_mono, , ] %>% tf$reduce_sum(axis = c(1L, 3L)) %>% as.matrix()

	}


	gr

} # predict.vae_baseline

