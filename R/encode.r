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

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		z <- matrix(NA, length(x), model@latent_dim)

		for (i in 1:n_batch){

			if (i %% 10 == 0)
				flog.info(sprintf('predicting | batch=%5.d/%5.d', i, n_batch))

			b <- starts[i]:ends[i]

			z[b, ]<- x[b]$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					metadata(x)$n_bins_per_window, metadata(x)$n_intervals,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				model@encoder() %>%
				as.matrix()

		}

		z
	}
) # encode


#' predict.vae_baseline
#'
setMethod(
	'encode',
	signature(
		model = 'vplot_knn_autoencoder_model',
		x = 'GRanges'
	),
	function(
		model,
		x,
		batch_size = 128   # v-plot per batch
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

			z[b, ]<- x[b]$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					metadata(x)$n_bins_per_window, metadata(x)$n_intervals,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				model@encoder() %>%
				as.matrix()

		}

		z
	}
) # encode


setMethod(
	'encode',
	signature(
		model = 'vplot_autoencoder_model',
		x = 'tensorflow.tensor'
	),
	function(
		model,
		x,
		batch_size = 128   # v-plot per batch
	){


		n <- dim(x)[1]
		starts <- seq(1, n, by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > n] <- n
		n_batch <- length(starts)

		z <- matrix(NA, n, model@latent_dim)

		for (i in 1:n_batch){

			if (i %% 10 == 0)
				flog.info(sprintf('predicting | batch=%5.d/%5.d', i, n_batch))

			b <- starts[i]:ends[i]

			z[b, ]<- x[b, , , ] %>%
				model@encoder() %>%
				as.matrix()
		}

		z
	}
) # encode


setMethod(
	'encode',
	signature(
		model = 'vplot_knn_autoencoder_model',
		x = 'dgCMatrix'
	),
	function(
		model,
		x,
		batch_size = 128,   # v-plot per batch
		n_bins_per_window,
		n_intervals
	){

		n <- nrow(x)
		starts <- seq(1, n, by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > n] <- n
		n_batch <- length(starts)

		z <- matrix(NA, n, model@latent_dim)

		for (i in 1:n_batch){

			if (i %% 10 == 0)
				flog.info(sprintf('predicting | batch=%5.d/%5.d', i, n_batch))

			b <- starts[i]:ends[i]

			z[b, ]<- x[b, ] %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					n_bins_per_window, n_intervals,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				model@encoder() %>%
				as.matrix()

		}

		z
	}
) # encode
