#' predict.vae_baseline
#'
setMethod(
	'encode',
	signature(
		model = 'vplot_autoencoder_model',
		x = 'GRanges'
	),
	function(
		model,
		x,
		batch_size = 128   # v-plot per batch
	){

		interval_per_batch <- floor(batch_size / metadata(x)$n_samples)

		starts <- seq(1, length(x), by = interval_per_batch)
		ends <- starts + interval_per_batch - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		z <- array(NA, dim = c(length(x), metadata(x)$n_samples, model@latent_dim))

		for (i in 1:n_batch){

			if (i %% 10 == 0)
				flog.info(sprintf('predicting | batch=%5.d/%5.d', i, n_batch))

			b <- starts[i]:ends[i]
			batch_size2 <- length(b) * metadata(x)$n_samples

			xi <- mcols(x[b])$counts %>%
				array_reshape(c(
					length(b),
					metadata(x)$n_bins_per_window,
					metadata(x)$n_intervals,
					metadata(x)$n_samples
				)) %>%
				array_permute(c(1, 4, 2, 3)) %>%
				array_reshape(c(
					batch_size2,
					metadata(x)$n_bins_per_window * metadata(x)$n_intervals
				)) %>%
				as_dgCMatrix() %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					batch_size2,
					metadata(x)$n_bins_per_window, metadata(x)$n_intervals,
					1L
					)) %>%
				tf$cast(tf$float32)

			zi <- xi %>%
				model@encoder()
				
			zi <- zi$mean() %>%
				as.matrix()

			dim(zi) <- c(length(b), metadata(x)$n_samples, model@latent_dim)
			z[b, , ] <- zi
		}

		z <- z %>% aperm(c(1, 3, 2))
		dim(z) <- c(length(x), metadata(x)$n_samples * model@latent_dim)
		z
	}
) # encode
