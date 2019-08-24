#' vae
#'

vae <- function(input_dim, feature_dim, latent_dim, window_size, num_samples){

	vplot_input <- layer_input(shape = c(feature_dim, input_dim, 1L))
	coverage_input <- layer_input(shape = c(window_size, 1L))

	z_vplot <- vplot_input %>% vplot_encoder_model()
	z_coverage <- coverage_input %>% coverage_encoder_model()

	z <- layer_concatenate(list(z_vplot, z_coverage)) %>%
	  layer_dense(units = params_size_multivariate_normal_tri_l(latent_dim)) %>%
	  layer_multivariate_normal_tri_l(event_size = latent_dim) %>%
		layer_kl_divergence_add_loss(
			distribution = tfd_independent(
				tfd_normal(loc = rep(0, latent_dim), scale = 1),
				reinterpreted_batch_ndims = 1
			),
		weight = 1)

	vplot_decoded <- z %>% vplot_decoder_model(input_dim = input_dim, feature_dim = feature_dim)
	vplot_prob <- vplot_decoded %>% 
		layer_flatten() %>%
		layer_independent_bernoulli(
			event_shape = c(feature_dim, input_dim, 1L),
			convert_to_tensor_fn = tfp$distributions$Bernoulli$logits,
			name = 'vplot_output'
		)

	coverage_decoded <- z %>% coverage_decoder_model(window_size = window_size)
	coverage_prob <- coverage_decoded %>%
		layer_flatten() %>%
		layer_independent_bernoulli(
			event_shape = c(window_size, 1L),
			convert_to_tensor_fn = tfp$distributions$Bernoulli$logits,
			name = 'coverage_output'
		)

	structure(list(
		vae = keras_model(
			inputs = list(vplot_input, coverage_input),
			outputs = list(vplot_prob, coverage_prob)
		),
		decoded = keras_model(
			inputs = list(vplot_input, coverage_input),
			outputs = list(vplot_decoded, coverage_decoded)
		),
		latent = keras_model(
			inputs = list(vplot_input, coverage_input),
			outputs = z
		),
		num_samples = num_samples,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		window_size = window_size
	), class = 'vae')

} # vae


#' save_model
#'
save_model <- function(x, dir){

	if (missing(dir))
		stop('dir must be specified')

	if (!file.exists(dir))
		dir.create(dir, recursive = TRUE)

	vae_file <- sprintf('%s/vae.h5', dir)
	flog.info(sprintf('writing %s', vae_file))
	x$vae$save_weights(vae_file)
	x$vae_file <- vae_file
	x$vae <- NULL

	latent_file <- sprintf('%s/latent.h5', dir)
	flog.info(sprintf('writing %s', latent_file))
	x$latent$save_weights(latent_file)
	x$latent_file <- latent_file
	x$latent <- NULL

	decoded_file <- sprintf('%s/decoded.h5', dir)
	flog.info(sprintf('writing %s', decoded_file))
	x$decoded$save_weights(decoded_file)
	x$decoded_file <- decoded_file
	x$decoded <- NULL

	model_file <- sprintf('%s/model.rds', dir)
	flog.info(sprintf('writing %s', model_file))
	saveRDS(x, model_file)

} # save_model


#' load_model
load_model <- function(dir){

	if (missing(dir))
		stop('dir must be specified')

	model_file <- sprintf('%s/model.rds', dir)

	if (!file.exists(model_file))
		stop(sprintf('%s does not exist', model_file))

	flog.info(sprintf('reading %s', model_file))
	x <- readRDS(model_file)

	model <- vae(
		input_dim = x$input_dim,
		feature_dim = x$feature_dim,
		latent_dim = x$latent_dim,
		window_size = x$window_size,
		num_samples = x$num_samples
	)

	flog.info(sprintf('reading %s', x$vae_file))
	model$vae$load_weights(x$vae_file)

	flog.info(sprintf('reading %s', x$latent_file))
	model$latent$load_weights(x$latent_file)

	flog.info(sprintf('reading %s', x$decoded_file))
	model$decoded$load_weights(x$decoded_file)

	model

} # load_model



#' vplot_decoder_model
#'
#' @param input_dim input dimension of the v-plot (i.e. the width of the genomic region)
#' @param feature_dim feature dimension of the v-plot	(i.e. the fragment size dimension)
#' @param filters0 the beginning filter dimension coming out of the latent layer
#' @param filters the filter sizes of each deconv layer
#' @param kernel_size the kernel size of each deconv layer.  The feature and input spaces shared the same kernel size. 
#'
vplot_decoder_model <- function(
	x,
	input_dim, 
	feature_dim, 
	filters0 = 64L, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 2L)
){

	input_dim0 <- input_dim / prod(input_strides)
  feature_dim0 <- feature_dim / prod(feature_strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	y <- x %>% 
		layer_dense(units = output_dim0, activation = 'relu') %>%
		layer_dropout(rate = 0.2) %>%
		layer_reshape(target_shape = c(feature_dim0, input_dim0, filters0)) %>%

		layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(feature_strides[1], input_strides[1]),
			padding = 'same',
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%

		layer_conv_2d_transpose( 
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(feature_strides[2], input_strides[2]),
			padding = 'same',
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%

		layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(feature_strides[3], input_strides[3]),
			padding = 'same',
			name = 'vplot_decoded'
		) 
	y
} # vplot_decoder_model


coverage_decoder_model <- function(
	x,
	window_size,
	filters0 = 16L, 
	filters = c(8L, 1L), 
	kernel_size = c(3L, 3L), 
	strides = c(4L, 1L))
{

	window_size0<- window_size / prod(strides)
	output_dim0 <- window_size0 * 1 * filters0

	y <- x %>%
		layer_dense(units = output_dim0, activation = 'relu') %>%
		layer_dropout(rate = 0.2) %>%
		layer_reshape(target_shape = c(window_size0, 1L, filters0)) %>%

		layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = shape(kernel_size[1], 1L),
			strides = shape(strides[1], 1L),
			padding = 'same',
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%

		layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = shape(kernel_size[2], 1L),
			strides = shape(strides[2], 1L),
			padding = 'same'
		) %>% 
		layer_reshape(
			target_shape = c(window_size, 1L),
			name = 'coverage_decoded'
		)
	y			
} # coverage_decoder_model


#' vplot_encoder_model
#'
vplot_encoder_model <- function(x, output_dim = 16L){

	y <- x %>%
		layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_flatten() %>%
		layer_dense(units = output_dim, activation = 'relu')

} # vplot_encoder_model


coverage_encoder_model <- function(x, output_dim = 8L){

	y <- x %>%
		layer_conv_1d(
			filters = 8L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_1d(
			filters = 8L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_1d(
			filters = 8L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_flatten() %>%
		layer_dense(units = output_dim, activation = 'relu')

} # coverage_encoder_model


