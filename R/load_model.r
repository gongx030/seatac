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
		slot_names <- slot_names[!slot_names %in% c('encoder', 'decoder')]
		param <- lapply(slot_names, function(x) slot(model, x))
		names(param) <- slot_names

		model <- do.call(build_model, c(name = class(model), x = NULL, param))

		# initialize the weights
		k_random_uniform(c(1L, model@n_bins_per_window, model@n_intervals, 1L)) %>%
			model@encoder() %>%
			model@decoder()

		encoder_file <- sprintf('%s/encoder.h5', dir)
		model@encoder$load_weights(encoder_file)

		decoder_file <- sprintf('%s/decoder.h5', dir)
		model@decoder$load_weights(decoder_file)

		model	

	}
)
