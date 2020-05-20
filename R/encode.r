setMethod(
	'encode',
	signature(
		model = 'vplot_autoencoder_model'
	),
	function(
		model,
		batch_size = 128   # v-plot per batch
	){

		starts <- seq(1, length(model@data), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(model@data)] <- length(model@data)
		n_batch <- length(starts)

		z <- matrix(NA, length(model@data), model@latent_dim)

		for (i in 1:n_batch){

			b <- starts[i]:ends[i]

			z[b, ]<- model@data[b]$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					metadata(model@data)$n_bins_per_window, metadata(model@data)$n_intervals,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				model@encoder() %>%
				as.matrix()

		}

		z
	}
) # encode

