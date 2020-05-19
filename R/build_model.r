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

		if (name == 'vplot_autoencoder'){

			latent_dim <- param[['latent_dim']]

			flog.info(sprintf('latent dim(latent_dim):%d', latent_dim))

			model <- new(
				'vplot_autoencoder_model',
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
		}else
			stop(sprintf('unknown name: %s', name))

		model
	}
)
