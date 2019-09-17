
#' vae
vae <- function(input_dim, feature_dim, latent_dim, num_samples){

	vplot_input <- layer_input(shape = c(feature_dim, input_dim, 1L))

	z <- vplot_input %>%
		vplot_encoder_model() %>%
		layer_dense(units = params_size_multivariate_normal_tri_l(latent_dim)) %>%
		layer_multivariate_normal_tri_l(event_size = latent_dim) %>%
		layer_kl_divergence_add_loss(
			distribution = tfd_independent(
			tfd_normal(loc = rep(0, latent_dim), scale = 1),
			reinterpreted_batch_ndims = 1
			)
		)

	vplot_decoded <- z %>% vplot_decoder_model(input_dim = input_dim, feature_dim = feature_dim)

	vplot_prob <- vplot_decoded %>% 
		layer_flatten() %>%
		layer_independent_bernoulli(
			event_shape = c(feature_dim, input_dim, 1L),
			convert_to_tensor_fn = tfp$distributions$Bernoulli$logits
		)

	structure(list(
		vae = keras_model(
			inputs = vplot_input,
			outputs = vplot_prob
		),
		decoded = keras_model(
			inputs = vplot_input,
			outputs = vplot_decoded
		),
		latent = keras_model(
			inputs = vplot_input,
			outputs = z
		),
		num_samples = num_samples,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim
	), class = 'vae')

} # vae

