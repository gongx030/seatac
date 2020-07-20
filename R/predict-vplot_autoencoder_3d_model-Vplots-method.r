#' vplot_parametric_vae_model
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_autoencoder_3d_model',
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

		latent <- array(NA, c(length(x),  model@n_blocks_per_window, model@latent_dim))
		predicted_counts <- Matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			b <- starts[i]:ends[i]

			xi <- x[b] %>%
				prepare_vplot(model) %>%
				extract_blocks_from_vplot(model@n_bins_per_block)

			w <- xi %>% tf$reduce_sum(c(2L, 3L), keepdims = TRUE)
			xi <- xi %>% tf$multiply( 1 / w)

			z <- xi %>%
				model@encoder()

			xi_pred <- z %>%
				model@decoder() %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			xi_pred <- aperm(xi_pred, c(1, 3, 2))
			dim(xi_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)
			xi_pred <- as(xi_pred, 'dgCMatrix')

			latent[b, , ] <- z %>%
				as.matrix()

			predicted_counts[b, ] <- xi_pred

		}

		mcols(x)$latent <- latent
		mcols(x)$predicted_counts <- predicted_counts

		x
	}
) # predict

