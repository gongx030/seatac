#' vplot_parametric_vae_model
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_vae_model',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 32L # v-plot per batch
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		latent <- array(NA, c(length(x),  model@n_blocks_per_window, model@latent_dim))
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			b <- starts[i]:ends[i]

			xi <- x[b] %>%
				prepare_vplot() %>%
				extract_blocks_from_vplot(model@n_bins_per_block) %>%
				tf$clip_by_value(clip_value_min = 0, clip_value_max = model@max_reads_per_pixel) %>%
				tf$reshape(c(length(b) * model@n_blocks_per_window, model@n_intervals, model@n_bins_per_block, 1L))

			posterior <- xi %>%
				model@encoder()

			zi <- posterior$mean() 

			xi_pred <- zi %>%
				model@decoder()

			xi_pred <- xi_pred$mean() %>%
				tf$reshape(c(length(b), model@n_blocks_per_window, model@n_intervals, model@n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			xi_pred <- aperm(xi_pred, c(1, 3, 2))
			dim(xi_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)

			zi <- zi %>%
				tf$reshape(c(length(b), model@n_blocks_per_window, model@latent_dim))

			latent[b, , ] <- zi %>%
				as.array()

			predicted_counts[b, ] <- xi_pred

		}

		mcols(x)$latent <- latent
		mcols(x)$predicted_counts <- predicted_counts
		class(x) <- 'VplotsFitted'
		x@model <- model

		x
	}
) # predict

