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

#' sequence_encoder_model
#'
sequence_encoder_model <- function(x, output_dim = 16L){

	y <- x %>%
		layer_conv_1d(
			filters = 32L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_1d(
			filters = 32L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_conv_1d(
			filters = 32L,
			kernel_size = 3L,
			strides = 2L,
			activation = 'relu'
		) %>%
		layer_batch_normalization() %>%
		layer_flatten() %>%
		layer_dense(units = output_dim, activation = 'relu')

} # sequence_encoder_model

#' sequence_encoder_model
#'
sequence_encoder_model2 <- function(x, output_dim = 16L){

	y <- x %>%
		bidirectional(layer_cudnn_gru(units = output_dim / 2)) %>%
		layer_dropout(0.2)
	y

} # sequence_encoder_model

