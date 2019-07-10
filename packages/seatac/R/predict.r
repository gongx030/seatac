predict.vae <- function(model, x, batch_size = 256){

	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples
	sample_dim <- window_dim * num_samples

	X <- mcols(x)$counts %>%
		as.matrix() %>%
		tf$cast(tf$float32) %>% 
		tf$reshape(shape(window_dim, num_samples, model$input_dim, model$feature_dim)) %>%
		tf$reshape(shape(window_dim * num_samples, model$input_dim, model$feature_dim)) %>%
		tf$expand_dims(axis = 3L)

	Z <- model$encoder(X)$loc

	if (model$prior == 'gmm'){
		P <- model$latent_prior_model(NULL)$components_distribution$log_prob(
			Z %>% tf$reshape(shape(sample_dim, 1, model$latent_dim))
		) %>% as.matrix()
	}else if (model$prior == 'hmm'){

		# See how to get the posterior of HMM
		browser()

		model$latent_prior_model(NULL)$posterior_marginals(
			Z %>% tf$reshape(shape(window_dim, num_samples, 1, model$latent_dim, 1))
		)

		browser()
	}

	t(matrix(max.col(P), num_samples, window_dim))

} # predict.vae
