#' predict
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_parametric_vae_v5_model',
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

		is_mono_nucleosome <- matrix(NA, length(x), x@n_bins_per_window)
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

			xi <- xi / tf$reduce_sum(xi, c(1L, 2L, 3L), keepdims = TRUE)

			posterior <- xi %>%
				model@encoder()

			z <- posterior$mean()

			v <- z %>%
				model@decoder()

			is_mono_nucleosome[b, ] <-  v %>% 
				as.matrix()

			latent[b, ] <- z %>% 
				as.matrix()

		}

		mcols(x)$is_mono_nucleosome <- is_mono_nucleosome
		mcols(x)$latent <- latent

		x
	}
) # predict

