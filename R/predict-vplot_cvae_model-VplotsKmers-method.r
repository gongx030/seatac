#'
setMethod(
	'predict',
	signature(
		model = 'vplot_cvae_model',
		x = 'VplotsKmers'
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

		block_size <- model@encoder$block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		latent <- array(NA, c(length(x),  n_blocks_per_window, model@encoder$latent_dim))
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting %s | batch=%4.d/%4.d', class(model), i, n_batch))

			b <- starts[i]:ends[i]

			inputs <- x[b] %>% prepare_blocks(model, min_reads = 0L)
			xi <- inputs$vplots
			ci <- inputs$kmers

			enc <- xi %>% model@encoder(ci)

			z <- enc$posterior$mean() 

			xi_pred <- z %>% model@decoder(enc$context)

			xi_pred <- xi_pred$mean() %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			xi_pred <- aperm(xi_pred, c(1, 3, 2))
			dim(xi_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)

			z <- z %>%
				tf$reshape(c(length(b), n_blocks_per_window, model@encoder$latent_dim))

			latent[b, , ] <- z %>%
				as.array()

			predicted_counts[b, ] <- xi_pred

		}

		mcols(x)$latent <- latent
		mcols(x)$predicted_counts <- predicted_counts

		x
	}
) # predict

