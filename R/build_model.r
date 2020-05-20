setMethod(
	'build_model',
	signature(
		x = 'GRanges',
		name = 'character'
	),
	function(
		x,
		name,
		...
	){

		param <- list(...)

		if (name == 'vplot_autoencoder_model'){

			latent_dim <- param[['latent_dim']]

			flog.info(sprintf('latent dim(latent_dim):%d', latent_dim))

			model <- new(
				name,
				encoder = encoder_model(
					latent_dim = latent_dim,
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L),
				),
				decoder = decoder_model(
					window_dim = metadata(x)$n_bins_per_window,
					interval_dim = metadata(x)$n_intervals,
					filters0 = 64,
					filters = c(32L, 32L, 1L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 2L),
				),
				n_samples = as.integer(metadata(x)$n_samples),
				window_dim = as.integer(metadata(x)$n_bins_per_window),
				interval_dim = as.integer(metadata(x)$n_intervals),
				latent_dim = as.integer(latent_dim)
			)

		}else if (name == 'vplot_autoencoder_cluster_model'){

			latent_dim <- param[['latent_dim']]
			num_clusters <- param[['num_clusters']]
			gamma <- param[['gamma']]
			sigma <- param[['sigma']]

			model <- new(
				name,
				encoder = encoder_model(
					latent_dim = latent_dim,
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L),
				),
				decoder = decoder_model(
					window_dim = metadata(x)$n_bins_per_window,
					interval_dim = metadata(x)$n_intervals,
					filters0 = 64,
					filters = c(32L, 32L, 1L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 2L),
				),
				n_samples = as.integer(metadata(x)$n_samples),
				window_dim = as.integer(metadata(x)$n_bins_per_window),
				interval_dim = as.integer(metadata(x)$n_intervals),
				latent_dim = as.integer(latent_dim),
				num_clusters = as.integer(num_clusters),
				gamma = gamma,
				sigma = sigma
			)

		}else if (name == 'vplot_autoencoder_rge_model'){

			latent_dim <- param[['latent_dim']]
			num_clusters <- param[['num_clusters']]
			gamma <- param[['gamma']]
			sigma <- param[['sigma']]
			lambda <- param[['lambda']]

			model <- new(
				name,
				encoder = encoder_model(
					latent_dim = latent_dim,
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L),
				),
				decoder = decoder_model(
					window_dim = metadata(x)$n_bins_per_window,
					interval_dim = metadata(x)$n_intervals,
					filters0 = 64,
					filters = c(32L, 32L, 1L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 2L),
				),
				n_samples = as.integer(metadata(x)$n_samples),
				window_dim = as.integer(metadata(x)$n_bins_per_window),
				interval_dim = as.integer(metadata(x)$n_intervals),
				latent_dim = as.integer(latent_dim),
				num_clusters = as.integer(num_clusters),
				gamma = gamma,
				sigma = sigma,
				lambda = lambda
			)
		}else
			stop(sprintf('unknown name: %s', name))

		model@data <- x

		model
	}
)
