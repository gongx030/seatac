#' predict.vae
#' 
predict.vae <- function(model, x, batch_size = 2^13){

	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples
	n_bins_per_window <- metadata(x)$window_size / metadata(x)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	# determine the batches
	bs <- seq(1, window_dim, by = batch_size)
	be <- bs + batch_size - 1
	be[be > window_dim] <- window_dim
	n_batch <- length(bs)

	latent <- NULL
	fitted_coverage  <- NULL
	fitted_counts <- NULL
	
	if (model$prior == 'gmm'){

		for (b in 1:n_batch){

			flog.info(sprintf('prediction | batch=%4.d/%4.d', b, n_batch))

			i <- bs[b]:be[b]
			window_dim2 <- length(i)

			X <- mcols(x)$counts[i, ] %>%
				as.matrix() %>%
				array_reshape(c(window_dim2, model$input_dim, model$feature_dim, 1L)) %>%
				tf$cast(tf$float32) # need to convert to tensors

			Y <- mcols(x)$coverage[i, , drop = FALSE] %>%
				array_reshape(c(window_dim2, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			label_pred <- model$classifier(list(vplot = X, coverage = Y))$mean() %>% as.matrix()

			browser()

			Z_x <- model$encoder$vplot(X)$loc
			Z_y <- model$encoder$coverage(Y)$loc
			Z <- tf$concat(list(Z_x, Z_y), axis = 1L)

			X <- model$decoder$vplot(Z)$mean() %>% 
				tf$squeeze() %>%
				as.array() %>%
				aperm(perm = c(1, 3, 2))	# change from row-major to column-major

			dim(X) <- c(window_dim2, model$input_dim * model$feature_dim)

			Y <- model$decoder$coverage(Z)$mean() %>% 
				as.array() 
	
			dim(Y) <- c(window_dim2, model$input_dim)

			latent <- rbind(latent, Z %>% as.matrix())
			fitted_coverage <- rbind(fitted_coverage, Y)
			fitted_counts <- rbind(fitted_counts, X)
		}
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))

	mcols(x)$latent <- latent
	mcols(x)$fitted_coverage <- fitted_coverage
	mcols(x)$fitted_counts <- fitted_counts
	x

} # predict.vae


