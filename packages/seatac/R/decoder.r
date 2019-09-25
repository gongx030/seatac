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
			padding = 'same'
		) 
	y
} # vplot_decoder_model
