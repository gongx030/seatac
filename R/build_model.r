setMethod(
	'build_model',
	signature(
		name = 'character',
		x = 'GRanges'
	),
	function(
		name,
		x,
		...
	){

		param <- list(...)

		if (name == 'vplot_autoencoder_model'){

			model <- build_model(
				name = name,
				n_bins_per_window = metadata(x)$n_bins_per_window,
				n_intervals = metadata(x)$n_intervals,
				n_samples = metadata(x)$n_samples,
				fragment_size_range = metadata(x)$fragment_size_range,
				fragment_size_interval = metadata(x)$fragment_size_interval,
				window_size = metadata(x)$window_size,
				bin_size = metadata(x)$bin_size,
				...
			)

		}else if (name == 'vplot_parametric_vae_model'){

			nfr_prob <- build_fragment_size_model(x, n = 1e4, bootstrap = 20)

			model <- build_model(
				name = name,
				n_bins_per_window = metadata(x)$n_bins_per_window,
				n_intervals = metadata(x)$n_intervals,
				n_samples = metadata(x)$n_samples,
				fragment_size_range = metadata(x)$fragment_size_range,
				fragment_size_interval = metadata(x)$fragment_size_interval,
				window_size = metadata(x)$window_size,
				bin_size = metadata(x)$bin_size,
				nfr_prob = nfr_prob,
				...
			)

			flog.info('building fragment size model')

		}else
			stop(sprintf('unknown name: %s', name))

		model
	}
)




setMethod(
	'build_model',
	signature(
		name = 'character',
		x = 'missing'
	),
	function(
		name,
		x,
		...
	){

		param <- list(...)

		if (name == 'vplot_autoencoder_model'){

			model <- new(
				name,
				encoder = encoder_model(
					latent_dim = as.integer(param$latent_dim),
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L),
				),
				decoder = decoder_model(
					window_dim = as.integer(param$n_bins_per_window),
					interval_dim = as.integer(param$n_intervals),
					filters0 = 64,
					filters = c(32L, 32L, 1L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 2L),
				),
				n_samples = as.integer(param$n_samples),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size)
			)
		}else if (name == 'vplot_parametric_vae_model'){

			param <- list(...)

			model <- new(
				name,
				encoder = vae_encoder_model(
					latent_dim = as.integer(param$latent_dim),
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L)
				),
				prior = prior_model(
					latent_dim = as.integer(param$latent_dim)
				),
				decoder = parametric_vae_decoder_model(
				),
				n_samples = as.integer(param$n_samples),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				nfr_prob = param$nfr_prob
			)

		}else if (name == 'vplot_autoencoder_cluster_v2_model'){

			model <- do.call(build_model, c(name = 'vplot_autoencoder_model', x = NULL, param))
			class(model) <- 'vplot_autoencoder_cluster_v2_model'
			model@num_clusters <- as.integer(param$num_clusters)
			model@sigma <- param$sigma
			model@gamma <- param$gamma
			model@membership <- NULL
			model@centers <- NULL

		}else
			stop(sprintf('unknown name: %s', name))

		model
	}
)

