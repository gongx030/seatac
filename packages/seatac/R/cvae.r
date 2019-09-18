#' cvae
#'
cvae <- function(input_dim, feature_dim, latent_dim, num_samples, window_size, sequence_dim){

	vplot_input <- layer_input(shape = c(feature_dim, input_dim, 1L))
	sequence_input <- layer_input(shape = window_size)

	h <- vplot_input %>%
		vplot_encoder_model(output_dim = sequence_dim)

	cond <- sequence_input %>%
		sequence_encoder_model(window_size = window_size, output_dim = sequence_dim)

	z <- layer_add(list(h, cond)) %>%
		layer_dense(units = params_size_multivariate_normal_tri_l(latent_dim)) %>%
		layer_multivariate_normal_tri_l(event_size = latent_dim) %>%
		layer_kl_divergence_add_loss(
			distribution = tfd_independent(
			tfd_normal(loc = rep(0, latent_dim), scale = 1),
			reinterpreted_batch_ndims = 1
			)
		)

	z_cond <- cond %>% layer_dense(latent_dim)

	vplot_decoded <- layer_add(list(z, z_cond)) %>% 
		vplot_decoder_model(input_dim = input_dim, feature_dim = feature_dim)


	vplot_prob <- vplot_decoded %>% 
		layer_flatten() %>%
		layer_independent_bernoulli(
			event_shape = c(feature_dim, input_dim, 1L),
			convert_to_tensor_fn = tfp$distributions$Bernoulli$logits
		)

	structure(list(
		vae = keras_model(
			inputs = list(vplot_input, sequence_input),
			outputs = vplot_prob
		),
		latent = keras_model(
			inputs = list(vplot_input, sequence_input),
			outputs = list(z, cond, z_cond)
		),
		num_samples = num_samples,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		sequence_dim = sequence_dim
	), class = 'cvae')

} # vae

