#' predict.gmm_cvae
#' 
predict.gmm_cvae <- function(model, gr, batch_size = 512){

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
	h <- matrix(NA, window_dim, model$latent_dim)	# h
	nfr <- matrix(NA, window_dim, n_bins_per_window)
	mono_nucleosome <- matrix(NA, window_dim, n_bins_per_window)

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		sequence <- unlist(strsplit(as.character(mcols(gr[i])$sequence), '')) %>%
			factor(c('N', 'A', 'C', 'G', 'T')) %>%
			as.numeric() %>%
			matrix(nrow = length(i), ncol = window_size, byrow = TRUE)

		sequence <- sequence - 1

		g <- model$encoder(list(x, sequence))
		z[i, ] <- g[[1]]$mean() %>% as.matrix()
		h[i, ] <- g[[2]] %>% as.matrix()

#		zh <- model$decoder(cbind(z[i, ], h[i, ]))$mean()
		zh <- model$decoder(z[i, ] + h[i, ])$mean()

		nfr[i, ] <- zh[, which(metadata(gr)$nfr), , ] %>% 
			tf$reduce_mean(axis = 1L) %>%
			tf$squeeze() %>%
			as.matrix()

		mono_nucleosome[i, ] <- zh[, which(metadata(gr)$mono_nucleosome), , ] %>% 
			tf$reduce_mean(axis = 1L) %>%
			tf$squeeze() %>%
			as.matrix()
	}

	mcols(gr)$latent <- z
	mcols(gr)$h <- h 
	mcols(gr)$nfr <- nfr 
	mcols(gr)$mono_nucleosome <- mono_nucleosome

	gr

} # predict.gmm_cvae
