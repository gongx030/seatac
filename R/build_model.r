#' build_model
setMethod(
	'build_model',
	signature(
		name = 'character',
		x = 'VplotsList'
	),
	function(
		name,
		x,
		...
	){

		x <- unlist(x)	# to Vplots
		build_model(name, x, ...)

	}
)


#' build_model
setMethod(
	'build_model',
	signature(
		name = 'character',
		x = 'Vplots'
	),
	function(
		name,
		x,
		...
	){

		build_model(
			name = name,
			n_bins_per_window = x@n_bins_per_window,
			n_intervals = x@n_intervals,
			fragment_size_range = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			window_size = x@window_size,
			bin_size = x@bin_size,
			...
		)
	}
)



#' build_model
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
				decoder_window = parametric_vae_decoder_model(
					output_dim = as.integer(param$n_bins_per_window),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(2L, 2L, 1L)
				),
				decoder_interval = parametric_vae_decoder_model(
					output_dim = as.integer(param$n_intervals),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(2L, 2L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				sigma0 = 1
			)

		}else if (name == 'vplot_parametric_vae_v2_model'){

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
					output_dim = as.integer(param$n_bins_per_window),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(2L, 2L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				sigma0 = 1
			)

		}else if (name == 'vplot_parametric_vae_v3_model'){

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
					output_dim = as.integer(param$n_intervals),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(2L, 2L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				sigma0 = 1,
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
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

