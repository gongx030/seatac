#' decoder_model_vae_baseline_conv
#'
#' @param window_dim window dimension of the v-plot (i.e. the width of the genomic region)
#' @param interval_dim interval dimension of the v-plot	(i.e. the fragment size dimension)
#' @param filters0 the beginning filter dimension coming out of the latent layer
#' @param filters the filter sizes of each deconv layer
#' @param kernel_size the kernel size of each deconv layer.  The interval and window spaces shared the same kernel size. 
#'
decoder_model <- function(
	window_dim, 
	interval_dim, 
	filters0 = 64, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	window_strides = c(2L, 2L, 2L), 
	interval_strides = c(2L, 2L, 2L), 
	name = NULL
){

	window_dim0 <- window_dim / prod(window_strides)
	interval_dim0 <- interval_dim / prod(interval_strides)
	output_dim0 <- window_dim0 * interval_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(
			units = output_dim0, 
			kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			activity_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			activation = 'relu',
		)

		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(window_dim0, interval_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(window_strides[1], interval_strides[1]),
			padding = 'same',
			kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			activity_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(window_strides[2], interval_strides[2]),
			padding = 'same',
			kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			activity_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(window_strides[3], interval_strides[3]),
			kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			bias_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
			padding = 'same'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3()

			tfd_independent(
				tfd_poisson(log_rate = y), 
				reinterpreted_batch_ndims = 3L
			)
		}
	})
} # decoder_model_vae_baseline_conv


