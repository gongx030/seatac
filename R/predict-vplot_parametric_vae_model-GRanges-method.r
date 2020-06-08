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
		batch_size = 128,   # v-plot per batch
		n = 20L
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		d_interval_mean <- matrix(NA, length(x), metadata(x)$n_intervals)
		d_interval_sd <- matrix(NA, length(x), metadata(x)$n_intervals)
		d_window_mean <- matrix(NA, length(x), metadata(x)$n_bins_per_window)
		d_window_sd <- matrix(NA, length(x), metadata(x)$n_bins_per_window)
		latent <- matrix(NA, length(x), model@latent_dim)

		for (i in 1:n_batch){

			b <- starts[i]:ends[i]

			xi <- x[b]$smoothed_counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					metadata(x)$n_intervals,
					metadata(x)$n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32)

			posterior <- xi %>%
				model@encoder()

			z <- posterior$sample(n)

			v_window <- z %>%
				tf$reshape(shape(length(b) * n, model@latent_dim)) %>%
				model@decoder_window()

			v_interval <- z %>%
				tf$reshape(shape(length(b) * n, model@latent_dim)) %>%
				model@decoder_interval()

			d_window_mean[b, ] <- v_window %>% 
				tf$reshape(shape(n, length(b), metadata(x)$n_bins_per_window)) %>%
				tf$reduce_mean(0L) %>%
				as.matrix()

			d_window_sd[b, ] <- v_window %>% 
				tf$reshape(shape(n, length(b), metadata(x)$n_bins_per_window)) %>%
				tf$math$reduce_std(0L) %>%
				as.matrix()

			d_interval_mean[b, ] <- v_interval %>% 
				tf$reshape(shape(n, length(b), metadata(x)$n_intervals)) %>%
				tf$reduce_mean(0L) %>%
				as.matrix()

			d_interval_sd[b, ] <- v_interval %>% 
				tf$reshape(shape(n, length(b), metadata(x)$n_intervals)) %>%
				tf$math$reduce_std(0L) %>%
				as.matrix()

			latent[b, ] <- posterior$mean() %>%
				as.matrix()

		}

		mcols(x)$distance_window_mean <- d_window_mean
		mcols(x)$distance_window_sd <- d_window_sd
		mcols(x)$distance_interval_mean <- d_interval_mean
		mcols(x)$distance_interval_sd <- d_interval_sd
		mcols(x)$latent <- latent

		x
	}
) # predict

