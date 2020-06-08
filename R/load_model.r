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
		model <- readRDS(model_file)

		slot_names <- slotNames(model)
		slot_names <- slot_names[!slot_names %in% c('encoder', 'decoder', 'prior')]
		param <- lapply(slot_names, function(x) slot(model, x))
		names(param) <- slot_names

		model <- do.call(build_model, c(name = class(model), x = NULL, param))

		# initialize the weights
		z <- k_random_uniform(c(1L, model@n_intervals, model@n_bins_per_window, 1L)) %>%
			model@encoder() 
		
		z$sample() %>%
			model@decoder_window()

		z$sample() %>%
			model@decoder_interval()

		encoder_file <- sprintf('%s/encoder.h5', dir)
		model@encoder$load_weights(encoder_file)
		flog.info(sprintf('reading %s', encoder_file))

		prior_file <- sprintf('%s/prior.h5', dir)
		model@prior$load_weights(prior_file)
		flog.info(sprintf('reading %s', prior_file))

		decoder_window_file <- sprintf('%s/decoder_window.h5', dir)
		model@decoder_window$load_weights(decoder_window_file)
		flog.info(sprintf('reading %s', decoder_window_file))

		decoder_interval_file <- sprintf('%s/decoder_interval.h5', dir)
		model@decoder_interval$load_weights(decoder_interval_file)
		flog.info(sprintf('reading %s', decoder_interval_file))

		model	

	}
)
