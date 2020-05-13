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
		batch_size = 128,   # v-plot per batch
		recurrent_steps = 3
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		z <- matrix(NA, length(x), model@latent_dim)

		for (i in 1:n_batch){

			if (i %% 10 == 0)
				flog.info(sprintf('predicting | batch=%5.d/%5.d', i, n_batch))

			b <- starts[i]:ends[i]

			xi <- mcols(x[b])$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					metadata(x)$n_bins_per_window, metadata(x)$n_intervals,
					1L
				)) %>%
				tf$cast(tf$float32)
			
			w <- (xi > 0) %>% # index for the non-zero terms
				tf$cast(tf$float32)

			for (iter in seq_len(recurrent_steps)){

				if (iter == 1)
					xi_input <- xi
				else
					xi_input <- xi * w + xi_pred$mean() * (1 - w)

				posterior <- xi_input %>% model@encoder()
				posterior_sample <- posterior$sample()
				xi_pred <- posterior_sample %>% model@decoder()
			}

			z[b, ] <- posterior$mean() %>% as.matrix()

		}

		z
	}
) # encode
