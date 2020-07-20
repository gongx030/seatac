#' vplot_parametric_vae_model
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_autoencoder_cluster_model',
		x = 'Vplots'
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

		latent <- matrix(NA, length(x), model@latent_dim)

		for (i in 1:n_batch){

			b <- starts[i]:ends[i]

			xi <- x[b]$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					x@n_intervals,
					x@n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				tf$nn$conv2d(model@gaussian_kernel, strides = c(1, 1, 1, 1), padding = 'SAME')

			xi_min <- tf$reduce_min(xi, c(1L, 2L), keepdims = TRUE)
			xi_max <- tf$reduce_max(xi, c(1L, 2L), keepdims = TRUE)
			xi <- (xi - xi_min) / (xi_max - xi_min)

			zi <- xi %>%
				model@encoder()

			latent[b, ] <- zi %>%
				as.matrix()

		}

		mcols(x)$latent <- latent

		x
	}
) # predict

