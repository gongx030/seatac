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

		latent <- array(NA, c(length(x),  model@n_blocks_per_window, model@latent_dim))
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting %s | batch=%4.d/%4.d', class(model), i, n_batch))

			b <- starts[i]:ends[i]

			inputs <- x[b] %>% prepare_blocks(model, min_reads = 0L)
			xi <- inputs$vplots
			gi <- inputs$kmers

			ci <- gi %>%
				model@embedder()

			posterior <- list(xi, ci) %>%
				model@encoder()

			zi <- posterior$mean() 

			xi_pred <- list(zi, ci) %>%
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
		class(x) <- 'VplotsKmersFitted'
		x@model <- model

		x
	}
) # predict

