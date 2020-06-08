setMethod(
	'save_model',
	signature(
		x = 'vplot_parametric_vae_model',
		dir = 'character'
	),
	function(
		x,
		dir,
		...
	){

		if (!file.exists(dir))
			dir.create(dir, recursive = TRUE)

		encoder_file <- sprintf('%s/encoder.h5', dir)
		flog.info(sprintf('writing %s', encoder_file))
		x@encoder$save_weights(encoder_file)

		prior_file <- sprintf('%s/prior.h5', dir)
		flog.info(sprintf('writing %s', prior_file))
		x@prior$save_weights(prior_file)

		decoder_window_file <- sprintf('%s/decoder_window.h5', dir)
		flog.info(sprintf('writing %s', decoder_window_file))
		x@decoder_window$save_weights(decoder_window_file)

		decoder_interval_file <- sprintf('%s/decoder_interval.h5', dir)
		flog.info(sprintf('writing %s', decoder_interval_file))
		x@decoder_interval$save_weights(decoder_interval_file)

		model_file <- sprintf('%s/model.rds', dir)
		flog.info(sprintf('writing %s', model_file))
		saveRDS(x, model_file)

	}
)

