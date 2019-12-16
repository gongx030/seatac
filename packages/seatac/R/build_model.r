build_model <- function(type, ...){

	if (type == 'cvae'){
		model <- build_model_cvae(...)
	}else
		stop(sprintf('unknown model type: %s', type))

	model

} # build_model


#' cvae
#'
build_model_cvae <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_cvae(latent_dim, window_size),
		decoder = decoder_model_cvae(input_dim, feature_dim),
		latent_prior_model = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('cvae'))

} # vae

