#' vplot_parametric_vae_model
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_parametric_vae_v2_model',
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

		fs_mean <- matrix(NA, length(x), metadata(x)$n_bins_per_window)
		fs_sd <- matrix(NA, length(x), metadata(x)$n_bins_per_window)
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

			v <- z %>%
				tf$reshape(shape(length(b) * n, model@latent_dim)) %>%
				model@decoder()

			fs_mean[b, ] <-  v %>% 
				tf$reshape(shape(n, length(b), metadata(x)$n_bins_per_window)) %>%
				tf$reduce_mean(0L) %>%
				as.matrix()

			fs_sd[b, ] <- v %>% 
				tf$reshape(shape(n, length(b), metadata(x)$n_bins_per_window)) %>%
				tf$math$reduce_std(0L) %>%
				as.matrix()

			latent[b, ] <- posterior$mean() %>%
				as.matrix()

		}

		mcols(x)$fragment_size_mean <- fs_mean
		mcols(x)$fragment_size_sd <- fs_sd
		mcols(x)$latent <- latent

		x
	}
) # predict

