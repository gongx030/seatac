#' vplot_parametric_vae_model
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_parametric_vae_model',
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

		latent <- matrix(NA, length(x), model@latent_dim)
		distance_mean <- matrix(NA, length(x), x@n_intervals)
		position_mean <- matrix(NA, length(x), x@n_bins_per_window)

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

			xi <- xi / tf$reduce_sum(xi, c(1L, 2L, 3L), keepdims = TRUE)

			posterior <- xi %>%
				model@encoder()

			z <- posterior$mean() 

			v_window <- z %>%   # batch size ~ n_bins_per_window
				model@decoder_window()

			v_interval <- z %>%   # batch size ~ n_intervals
				model@decoder_interval()

			latent[b, ] <- z %>%
				as.matrix()

			distance_mean[b, ] <- v_interval %>%
				as.matrix()

			position_mean[b, ] <- v_window %>%
				as.matrix()

		}

		mcols(x)$latent <- latent
		mcols(x)$distance_mean <- distance_mean
		mcols(x)$position_mean <- position_mean

		x
	}
) # predict

