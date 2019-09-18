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
sequence_encoder_model <- function(x, window_size, output_dim = 16L){

	y <- x %>%
    layer_embedding(
			input_dim = 4L, 
			output_dim = 4L, 
			input_length = window_size, 
			weights = list(diag(4L)),
			trainable = FALSE
		) %>%
		layer_conv_1d(
			filters = 320L,
			kernel_size = 26L,
			strides = 1L,
			activation = 'relu'
		) %>%
		layer_max_pooling_1d(pool_size = 13L, strides = 13L) %>%
		layer_dropout(0.2) %>%
		bidirectional(layer_cudnn_gru(units = output_dim / 2)) %>%
		layer_dropout(0.5)
	y

} # sequence_encoder_model

