setMethod(
	'load_model',
	signature(
		dir = 'character'
	),
	function(
		dir,
		...
	){


		if (!file.exists(dir))
			stop(sprintf('%s does not exist dir'))

		model_file <- sprintf('%s/model.rds', dir)

		if (!file.exists(model_file))
			stop(sprintf('%s does not exist', model_file))

		flog.info(sprintf('reading %s', model_file))
		x <- readRDS(model_file)

		if (class(x) == 'vplot_autoencoder_model'){

			model <- x@data %>% build_model(name = 'vplot_autoencoder_model', latent_dim = x@latent_dim)

			# initialize the weights
			k_random_uniform(c(1L, model@window_dim, model@interval_dim, 1L)) %>%
				model@encoder() %>%
				model@decoder()

			encoder_file <- sprintf('%s/encoder.h5', dir)
			model@encoder$load_weights(encoder_file)

			decoder_file <- sprintf('%s/decoder.h5', dir)
			model@decoder$load_weights(decoder_file)

		}else
			stop(sprintf('unknown model class: %s', class(model)))

		model	

	}
)
