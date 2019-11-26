
#' vae
vae <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = gmm_cvae_encoder_model(latent_dim, sequence_dim, window_size),
		decoder = gmm_cvae_decoder_model(input_dim, feature_dim),
		latent_prior_model = gmm_prior_model(latent_dim, n_components),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		n_components = n_components,
		num_samples = num_samples,
		window_size = window_size,
		sequence_dim = sequence_dim
	), class = c('vae'))

} # vae

