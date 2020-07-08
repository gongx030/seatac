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
			activation = 'relu',
		)

		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(interval_dim0, window_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(interval_strides[1], window_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(interval_strides[2], window_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(interval_strides[3], window_strides[3]),
			padding = 'same',
			activation = 'relu'
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
			y

		}
	})
} # decoder_model_vae_baseline_conv


parametric_vae_decoder_model <- function(
	output_dim, 
	filters0 = 8L, 
	filters = c(8L, 8L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	interval_strides = c(2L, 2L, 2L), 
	name = NULL
){

	output_dim0 <- output_dim / prod(interval_strides)
  input_dim0 <- output_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(
			units = input_dim0,
			activation = 'relu'
		)

		self$dropout_1 <- layer_dropout(rate = 0.2)

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(1L, interval_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(1L, interval_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(1L, interval_strides[3]),
			padding = 'same',
			activation = 'relu'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
			  self$dropout_1() %>%
				layer_reshape(target_shape = c(1L, output_dim0, filters0)) %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3() %>%
				layer_reshape(target_shape = c(output_dim))

			y
		}
	})
}


parametric_vae_mixture_decoder_model <- function(
	output_dim, 
	filters0 = 8L, 
	filters = c(8L, 8L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	interval_strides = c(2L, 2L, 2L), 
	name = NULL
){

	output_dim0 <- output_dim / prod(interval_strides)
  input_dim0 <- output_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(
			units = input_dim0,
			activation = 'relu'
		)

		self$dropout_1 <- layer_dropout(rate = 0.2)

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(1L, interval_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(1L, interval_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(1L, interval_strides[3]),
			padding = 'same',
			activation = 'relu'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
			  self$dropout_1() %>%
				layer_reshape(target_shape = c(1L, output_dim0, filters0)) %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3() %>%
				layer_reshape(target_shape = c(output_dim))

			y
		}
	})
}

mixture_predictor_model <- function(
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(
			units = 1L,
			activation = 'sigmoid'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1()
			y
		}
	})
}

window_decoder_model <- function(
	output_dim, 
	filters0 = 8L, 
	filters = c(8L, 8L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	window_strides = c(2L, 2L, 2L), 
	name = NULL
){

	output_dim0 <- output_dim / prod(window_strides)
  input_dim0 <- output_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(
			units = input_dim0,
			activation = 'relu'
		)

		self$dropout_1 <- layer_dropout(rate = 0.2)

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(1L, window_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(1L, window_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(1L, window_strides[3]),
			padding = 'same',
			activation = 'sigmoid'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
			  self$dropout_1() %>%
				layer_reshape(target_shape = c(1L, output_dim0, filters0)) %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3() %>%
				layer_reshape(target_shape = c(output_dim))

			y
		}
	})
}

