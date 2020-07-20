#' predict
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_parametric_vae_v4_model',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 128,   # v-plot per batch
		n = 20L
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		distance_mean <- matrix(NA, length(x), x@n_intervals)
		latent <- matrix(NA, length(x), model@latent_dim)
		theta <- rep(NA, length(x))

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

			v <- z %>%
				model@decoder()

			theta[b] <- z %>%
				model@mixture() %>%
				as.matrix()

			distance_mean[b, ] <-  v %>% 
				as.matrix()

			latent[b, ] <- z %>% 
				as.matrix()

		}

		mcols(x)$distance_mean <- distance_mean
		mcols(x)$latent <- latent
		mcols(x)$theta <- theta

		x
	}
) # predict

