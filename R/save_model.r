#' save_model
#'
save_model <- function(x, dir){

	if (missing(dir))
		stop('dir must be specified')

	if (!file.exists(dir))
		dir.create(dir, recursive = TRUE)

	UseMethod('save_model', x)
}

#'
save_model.vae_baseline <- function(x, dir){

	encoder_file <- sprintf('%s/encoder.h5', dir)
	flog.info(sprintf('writing %s', encoder_file))
	x$encoder$save_weights(encoder_file)
	x$encoder_weight_file <- encoder_file

	decoder_file <- sprintf('%s/decoder.h5', dir)
	flog.info(sprintf('writing %s', decoder_file))
	x$decoder$save_weights(decoder_file)
	x$decoder_weight_file <- decoder_file

	prior_file <- sprintf('%s/prior.h5', dir)
	flog.info(sprintf('writing %s', prior_file))
	x$prior$save_weights(prior_file)
	x$prior_weight_file <- prior_file

	x$encoder <- NULL
	x$decoder <- NULL
	x$prior <- NULL

	model_file <- sprintf('%s/model.rds', dir)
	flog.info(sprintf('writing %s', model_file))
	saveRDS(x, model_file)

} # save_model.cvae
