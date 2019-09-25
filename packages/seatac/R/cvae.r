#' cvae
#'
cvae <- function(input_dim, feature_dim, latent_dim, num_samples, window_size, sequence_dim){

	filter_size <- 300L
	kernel_size <- 26L
	conv_strides <- 1L
	pool_size <- 13L
	pool_strides <- 13L

	vplot_input <- layer_input(shape = c(feature_dim, input_dim, 1L))
	sequence_input <- layer_input(shape = window_size)

	motifs <- sequence_input %>%
		layer_embedding(
			input_dim = 4L,
			output_dim = 4L,
			input_length = window_size,
			weights = list(diag(4L)),
			trainable = FALSE
		) %>% 
		layer_conv_1d(
			filters = filter_size,
			kernel_size = kernel_size,
			strides = conv_strides,
			activation = 'relu',
			name = 'sequence_conv',
		) 
		
	h_sequence <- motifs %>%
		layer_max_pooling_1d(pool_size = pool_size, strides = pool_strides) %>%
		layer_dropout(0.2) %>%
		bidirectional(layer_cudnn_gru(units = sequence_dim / 2)) %>%
		layer_dropout(0.5)
	
	h_vplot <- vplot_input %>%
		vplot_encoder_model(output_dim = sequence_dim)

	z <- layer_add(list(h_vplot, h_sequence)) %>%
		layer_dense(units = params_size_multivariate_normal_tri_l(latent_dim)) %>%
		layer_multivariate_normal_tri_l(event_size = latent_dim) %>%
		layer_kl_divergence_add_loss(
			distribution = tfd_independent(
			tfd_normal(loc = rep(0, latent_dim), scale = 1),
			reinterpreted_batch_ndims = 1
			)
		)

	z_sequence <- h_sequence %>%
		layer_dense(latent_dim)

	vplot_decoded <- layer_add(list(z, z_sequence)) %>% 
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
			outputs = z
		),
		sequence_motifs = keras_model(
			inputs = sequence_input,
			outputs = motifs
		),
		num_samples = num_samples,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		sequence_dim = sequence_dim,
		filter_size = filter_size,
		kernel_size = kernel_size,
		conv_strides = conv_strides,
		pool_size = pool_size,
		pool_strides = pool_strides
	), class = 'cvae')

} # vae

