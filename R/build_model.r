#' build_model
#'
build_model <- function(type = 'vae_baseline', input_dim, feature_dim, latent_dim, num_samples, window_size, ...){

	if(type == 'vae_baseline'){

		structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_baseline_conv(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'vplot',
			output = 'vplot',
			...
		), class = type)

	}else
		stop(sprintf('unknown model type: %s', type))

} # build_model

