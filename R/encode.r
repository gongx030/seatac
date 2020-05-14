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

		X <- x$counts
		X[X > 0] <- 1

		for (i in 1:n_batch){

			if (i %% 10 == 0)
				flog.info(sprintf('predicting | batch=%5.d/%5.d', i, n_batch))

			b <- starts[i]:ends[i]

			xi <- X[b, ] %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					metadata(x)$n_bins_per_window, metadata(x)$n_intervals,
					1L
				)) %>%
				tf$cast(tf$float32)
			
			for (iter in seq_len(recurrent_steps)){

				if (iter == 1)
					xi_input <- xi
				else
					xi_input <- xi * wi + xi_pred * (1 - wi)

				wi <- (xi_input > 0) %>%  # index for the non-zero terms
					tf$cast(tf$float32)

				zi <- xi_input %>% model@encoder()
				xi_pred <- zi %>% model@decoder()

			}

			if (i == 1)
				xi_pred %>% tf$reduce_sum(0L) %>% tf$squeeze() %>% as.matrix() %>% image()

			z[b, ] <- zi %>% as.matrix()

		}

		z
	}
) # encode
