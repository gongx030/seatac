
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
save_model.gmm_cvae <- function(x, dir){

	encoder_file <- sprintf('%s/encoder.h5', dir)
	flog.info(sprintf('writing %s', encoder_file))
	x$encoder$save_weights(encoder_file)
	x$encoder_weight_file <- encoder_file

	decoder_file <- sprintf('%s/decoder.h5', dir)
	flog.info(sprintf('writing %s', decoder_file))
	x$decoder$save_weights(decoder_file)
	x$decoder_weight_file <- decoder_file

	latent_prior_file <- sprintf('%s/latent_prior_model.h5', dir)
	flog.info(sprintf('writing %s', latent_prior_file))
	x$latent_prior_model$save_weights(latent_prior_file)
	x$latent_prior_weight_file <- latent_prior_file

	x$encoder <- NULL
	x$decoder <- NULL
	x$latent_prior_model <- NULL

	model_file <- sprintf('%s/model.rds', dir)
	flog.info(sprintf('writing %s', model_file))
	saveRDS(x, model_file)

} # save_model.gmm_cvae


#' load_model
load_model <- function(dir){

	if (missing(dir))
		stop('dir must be specified')

	model_file <- sprintf('%s/model.rds', dir)

	if (!file.exists(model_file))
		stop(sprintf('%s does not exist', model_file))

	flog.info(sprintf('reading %s', model_file))
	x <- readRDS(model_file)
	model <- initialize_gmm_cvae_model(x)

	for (s in names(x))
		if (is.null(model[[s]]))
			model[[s]] <- x[[s]]
	model

} # load_model


#' initialize_gmm_cvae_model
#'
initialize_gmm_cvae_model <- function(x){

	model <- gmm_cvae(
		input_dim = x$input_dim,
		feature_dim = x$feature_dim,
		latent_dim = x$latent_dim,
		n_components = x$n_components,
		num_samples = x$num_samples,
		window_size = x$window_size,
		sequence_dim = x$sequence_dim
	)

	vplot <- array(0, dim = c(1, model$feature_dim, model$input_dim, 1))  %>% tf$cast(tf$float32)
	sequence <- matrix(0, nrow = 1, ncol = model$window_size) %>% tf$cast(tf$float32)
	list(vplot, sequence) %>% model$encoder()
	model$encoder$load_weights(x$encoder_weight_file)

	z <- matrix(0, nrow = 1, ncol = model$latent_dim) %>% tf$cast(tf$float32)
	z %>% model$decoder()
	model$decoder$load_weights(x$decoder_weight_file)

	model$latent_prior_model$load_weights(x$latent_prior_weight_file)
	model

} # initialize_gmm_cvae_model

