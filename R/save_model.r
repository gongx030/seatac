setMethod(
	'save_model',
	signature(
		x = 'vplot_vae_model',
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

	}
)


setMethod(
	'save_model',
	signature(
		x = 'vplot_cvae_model',
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

		decoder_file <- sprintf('%s/decoder.h5', dir)
		flog.info(sprintf('writing %s', decoder_file))
		x@decoder$save_weights(decoder_file)

		embedder_file <- sprintf('%s/embedder.h5', dir)
		flog.info(sprintf('writing %s', embedder_file))
		x@embedder$save_weights(embedder_file)

		model_file <- sprintf('%s/model.rds', dir)
		flog.info(sprintf('writing %s', model_file))
		saveRDS(x, model_file)

	}
)
