#' gmm_cvae
#'
gmm_cvae <- function(input_dim, feature_dim, latent_dim, n_components, num_samples, window_size, sequence_dim, is_nfr, is_mono_nucleosome){

	flog.info(sprintf('num components(n_components):%d', n_components))
	flog.info(sprintf('sequence dimension(sequence_dim):%d', sequence_dim))

	structure(list(
		encoder = gmm_cvae_encoder_model(latent_dim, sequence_dim, window_size),
		decoder = gmm_cvae_decoder_model(input_dim, feature_dim),
		latent_prior_model = gmm_prior_model(latent_dim, n_components),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		n_components = n_components,
		num_samples = num_samples
	), class = c('gmm_cvae'))

} # vae


