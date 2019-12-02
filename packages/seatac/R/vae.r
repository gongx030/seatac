
#' vae
vae <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model(latent_dim, window_size),
		decoder = decoder_model(input_dim, feature_dim),
		latent_prior_model = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae'))

} # vae

