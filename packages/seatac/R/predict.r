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


encode <- function(model, x, window_size = 400, step_size = 200, batch_size = 2^12){

	expand <- metadata(x)$window_size
	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples
	bin_size <- metadata(x)$bin_size
	n_bins_per_window <- window_size / bin_size

	if (is.na(step_size))
		step_size <- window_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('input peak width: %d', expand))

	starts <- seq(1, expand / bin_size - n_bins_per_window + 1, by = step_size / bin_size)
	ends <- starts + n_bins_per_window - 1
	n_block <- length(starts)

	batch_size2 <- floor(batch_size / n_block)

	bs <- seq(1, window_dim, by = batch_size2)
	be <- bs + batch_size2 - 1
	be[be > window_dim] <- window_dim
	n_batch <- length(bs)

	gr <- NULL
									
	if (model$prior == 'gmm'){

		for (b in 1:n_batch){

			flog.info(sprintf('encoding | batch=%4.d/%4.d', b, n_batch))

			i <- bs[b]:be[b]
			xi <- makeData(x[i], window_size = window_size, step_size = step_size, min_reads_per_window = 0, min_reads_coverage = 0)
			window_dim2 <- length(xi)

			X <- mcols(xi)$counts %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>% 
				tf$reshape(shape(window_dim2, model$input_dim, model$feature_dim)) %>%
				tf$expand_dims(axis = 3L)

			Y <- mcols(xi)$coverage %>% 
				tf$cast(tf$float32) %>%
				tf$expand_dims(axis = 2L)

			Z_x <- model$encoder$vplot(X)$loc
			Z_y <- model$encoder$coverage(Y)$loc
			Z <- tf$concat(list(Z_x, Z_y), axis = 1L)

			Z <- Z %>% as.matrix()

			mcols(xi)$window_id <- i
			mcols(xi)$latent <- Z

			if (is.null(gr))
				gr <- xi
			else
				gr <- c(gr, xi)
		}
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))

	metadata(gr) <- metadata(x)
	metadata(gr)$window_size <- window_size
	metadata(gr)$n_bins_per_window <- n_bins_per_window
	gr	

} # encode


decode <- function(model, x, batch_size = 2^12){

	window_dim <- length(x)

	bs <- seq(1, window_dim, by = batch_size)
	be <- bs + batch_size - 1
	be[be > window_dim] <- window_dim
	n_batch <- length(bs)

	mcols(x)$fitted_counts <- matrix(0, length(x), model$input_dim * model$feature_dim)
	mcols(x)$fitted_coverage <- matrix(0, length(x), model$input_dim)

	if (model$prior == 'gmm'){

		for (b in 1:n_batch){

			flog.info(sprintf('decoding | batch=%4.d/%4.d', b, n_batch))
			i <- bs[b]:be[b]
			window_dim2 <- length(i)

			Z <- mcols(x[i])$latent

			mcols(x)$fitted_counts[i, ] <- model$decoder$vplot(Z)$mean() %>% 
				tf$squeeze() %>%
				tf$reshape(shape(window_dim2, model$input_dim * model$feature_dim)) %>%
				as.matrix() 

			mcols(x)$fitted_coverage[i, ] <- model$decoder$coverage(Z)$mean() %>% 
				tf$squeeze() %>%
				as.matrix() 
			
		}
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))
	x

} # decode


