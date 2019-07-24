predict.vae <- function(model, x, batch_size = 2^12){

	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples
	bin_size <- metadata(x)$bin_size

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)


	if (model$prior == 'gmm'){

		gr <- NULL

		for (b in 1:n_batch){

			flog.info(sprintf('prediction | batch=%4.d/%4.d', b, n_batch))

			i <- starts[b]:ends[b]
			xi <- x[i]
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

			mcols(xi)$fitted_counts <- model$decoder$vplot(Z)$mean() %>% 
				tf$squeeze() %>%
				tf$reshape(shape(window_dim2, model$input_dim * model$feature_dim)) %>%
				as.matrix() 

			mcols(xi)$fitted_coverage <- model$decoder$coverage(Z)$mean() %>% 
				tf$squeeze() %>%
				as.matrix() 

			if (is.null(gr))
				gr <- xi
			else
				gr <- c(gr, xi)
		}
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))

	gr
} # predict.vae
