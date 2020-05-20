setMethod(
	'save_model',
	signature(
		x = 'vplot_autoencoder_model',
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

		decoder_file <- sprintf('%s/decoder.h5', dir)
		flog.info(sprintf('writing %s', decoder_file))
		x@decoder$save_weights(decoder_file)

		model_file <- sprintf('%s/model.rds', dir)
		flog.info(sprintf('writing %s', model_file))
		saveRDS(x, model_file)

	}
)

