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
				bin_size = as.integer(param$bin_size),
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
			)
		}else if (name == 'vplot_autoencoder_v2_model'){

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
				bin_size = as.integer(param$bin_size),
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
			)

		}else if (name == 'vplot_vae_model'){

			model <- new(
				name,
				encoder = vae_encoder_model(
					latent_dim = as.integer(param$latent_dim),
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L),
				),
				prior = prior_model(
					latent_dim = as.integer(param$latent_dim)
				),
				decoder = decoder_model(
					window_dim = as.integer(param$n_bins_per_window),
					interval_dim = as.integer(param$n_intervals),
					filters0 = 64L,
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
				bin_size = as.integer(param$bin_size),
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
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
					interval_strides = c(1L, 1L, 1L)
				),
				decoder_interval = parametric_vae_decoder_model(
					output_dim = as.integer(param$n_intervals),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(1L, 1L, 1L)
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

			if (is.null(param$sigma0))
				param$sigma0 <- 1

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
					interval_strides = c(1L, 1L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				sigma0 = param$sigma0,
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
			)

		}else if (name == 'vplot_parametric_vae_v3_vampprior_model'){

			param <- list(...)

			if (is.null(param$sigma0))
				param$sigma0 <- 1

			if (is.null(param$num_pseudo_inputs))
				param$num_pseudo_inputs <- 2L

			model <- new(
				name,
				encoder = vae_encoder_model(
					latent_dim = as.integer(param$latent_dim),
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L)
				),
				prior = vampprior_prior_model(
					latent_dim = as.integer(param$latent_dim),
					num_pseudo_inputs = param$num_pseudo_inputs,
					window_dim = as.integer(param$n_bins_per_window),
					interval_dim = as.integer(param$n_intervals)
				),
				decoder = parametric_vae_decoder_model(
					output_dim = as.integer(param$n_intervals),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(1L, 1L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				sigma0 = param$sigma0,
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
			)

		}else if (name == 'vplot_parametric_vae_v3_gmm_model'){

			# GMM prior

			param <- list(...)

			if (is.null(param$sigma0))
				param$sigma0 <- 1

			if (is.null(param$n_components))
				param$n_components <- 5L

			model <- new(
				name,
				encoder = vae_encoder_model(
					latent_dim = as.integer(param$latent_dim),
					filters = c(32L, 32L, 32L),
					kernel_size = c(3L, 3L, 3L),
					window_strides = c(2L, 2L, 2L),
					interval_strides = c(2L, 2L, 1L)
				),
				prior = gmm_prior_model(
					latent_dim = as.integer(param$latent_dim),
					n_components = as.integer(param$n_components)
				),
				decoder = parametric_vae_decoder_model(
					output_dim = as.integer(param$n_intervals),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(1L, 1L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				sigma0 = param$sigma0,
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1)
			)

		}else if (name == 'vplot_parametric_vae_v4_model'){

			# Gaussian / uniform mixture model

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
				decoder = parametric_vae_mixture_decoder_model(
					output_dim = as.integer(param$n_intervals),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					interval_strides = c(2L, 2L, 1L)
				),
				mixture = mixture_predictor_model(),
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

		}else if (name == 'vplot_parametric_vae_v5_model'){

			# Assigning binary state to each 'bin' in the window

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
				decoder = window_decoder_model(
					output_dim = as.integer(param$n_bins_per_window),
					filters0 = 32L,
					filters = c(8L, 8L, 1L),
					kernel_size = c(6L, 6L, 6L),
					window_strides = c(1L, 1L, 1L)
				),
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size),
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1),
				mu = param$mu,
				sigma = param$sigma,
				shape = param$shape,
				scale = param$scale
			)

		}else if (name == 'vplot_autoencoder_cluster_model'){

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
				bin_size = as.integer(param$bin_size),
				gaussian_kernel = gaussian_kernel(size = 10, mu = 1, std = 1),
				num_clusters = as.integer(param$num_clusters),
				gamma = param$gamma,
				sigma = param$sigma

			)

		}else if (name == 'vplot_glove_model'){

			param <- list(...)

			model <- new(
				name,
				n_bins_per_window = as.integer(param$n_bins_per_window),
				n_intervals = as.integer(param$n_intervals),
				latent_dim = as.integer(param$latent_dim),
				fragment_size_range = as.integer(param$fragment_size_range),
				fragment_size_interval = as.integer(param$fragment_size_interval),
				window_size = as.integer(param$window_size),
				bin_size = as.integer(param$bin_size)
			)

		}else
			stop(sprintf('unknown name: %s', name))

		model
	}
)

