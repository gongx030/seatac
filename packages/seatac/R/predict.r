predict.vae <- function(model, x, batch_size = 2^14){

	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples
	sample_dim <- window_dim * num_samples

	X <- mcols(x)$counts %>%
		as.matrix() %>%
		tf$cast(tf$float32) %>% 
		tf$reshape(shape(window_dim, num_samples, model$input_dim, model$feature_dim)) %>%
		tf$reshape(shape(window_dim * num_samples, model$input_dim, model$feature_dim)) %>%
		tf$expand_dims(axis = 3L)

	starts <- seq(1, dim(X)[1], by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > dim(X)[1]] <- dim(X)[1]
	n_batch <- length(starts)

	if (model$prior == 'gmm'){
		P <- NULL
		for (b in 1:n_batch){
			flog.info(sprintf('prediction | batch=%4.d/%4.d', b, n_batch))
			i <- starts[b]:ends[b]
			window_dim_b <- length(i)
			Z <- model$encoder(X[i, , , , drop = FALSE])$loc %>%
				tf$reshape(shape(window_dim_b, 1, model$latent_dim))
			Pb <- model$latent_prior_model(NULL)$components_distribution$log_prob(Z) %>% as.matrix()
			P <- rbind(P, Pb)
		}
		C <- t(matrix(max.col(P), num_samples, window_dim))
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))

	C
} # predict.vae
