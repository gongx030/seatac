#' decoder_model_cvae
#'
#' @param input_dim input dimension of the v-plot (i.e. the width of the genomic region)
#' @param feature_dim feature dimension of the v-plot	(i.e. the fragment size dimension)
#' @param filters0 the beginning filter dimension coming out of the latent layer
#' @param filters the filter sizes of each deconv layer
#' @param kernel_size the kernel size of each deconv layer.  The feature and input spaces shared the same kernel size. 
#'
decoder_model_cvae <- function(
	input_dim, 
	feature_dim, 
	filters0 = 64, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 2L), 
	name = NULL
){

	input_dim0 <- input_dim / prod(input_strides)
	feature_dim0 <- feature_dim / prod(feature_strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(feature_dim0, input_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(input_strides[1], feature_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
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
					tfd_bernoulli(logits = y), 
					reinterpreted_batch_ndims = 3L
				)
		}
	})
} # decoder_model_cvae

#' decoder_model_vae
#'
#'
decoder_model_vae <- function(
	input_dim, 
	feature_dim, 
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = feature_dim)
		self$dense_2 <- layer_dense(units = input_dim)
		self$baseline <- tf$Variable(0)
			
		function(x, mask = NULL, training = TRUE){

			feature <- x %>% 
				self$dense_1() %>% 
				layer_dropout(rate = 0.2) %>%
				layer_reshape(target_shape = c(feature_dim, 1L, 1L))

			input <- x %>% 
				self$dense_2() %>% 
				layer_dropout(rate = 0.2) %>%
				layer_reshape(target_shape = c(1L, input_dim, 1L))

			y <- self$baseline + feature + input

			tfd_independent(
				tfd_bernoulli(logits = y),
				reinterpreted_batch_ndims = 3L
			)
		}
	})
} # decoder_model_vae


#' decoder_model_vae_baseline
#'
decoder_model_vae_baseline <- function(
	input_dim, 
	feature_dim, 
	filters0 = 64, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 2L), 
	name = NULL
){

	input_dim0 <- input_dim / prod(input_strides)
	feature_dim0 <- feature_dim / prod(feature_strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(feature_dim0, input_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(input_strides[1], feature_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
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
					tfd_bernoulli(logits = y), 
					reinterpreted_batch_ndims = 3L
				)
		}
	})
} # decoder_model_vae_baseline


#' decoder_model_vae_vplot_position_fragment_size
#'
decoder_model_vae_vplot_position_fragment_size <- function(
	input_dim, 
	feature_dim, 
	filters0 = 64, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 2L), 
	name = NULL
){

	input_dim0 <- input_dim / prod(input_strides)
	feature_dim0 <- feature_dim / prod(feature_strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(feature_dim0, input_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(input_strides[1], feature_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			padding = 'same'
		)

		self$dense_2 <- layer_dense(units = 8L, activation = 'relu')
		self$dense_3 <- layer_dense(units = feature_dim)

		self$dense_4 <- layer_dense(units = 8L, activation = 'relu')
		self$dense_5 <- layer_dense(units = feature_dim)

		function(x, mask = NULL){

			y1 <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3()

			y2 <- x %>%
				self$dense_2() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_3()

			y3 <- x %>%
				self$dense_4() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_5()

			list(
				vplot =	tfd_independent(
					tfd_bernoulli(logits = y1), 
					reinterpreted_batch_ndims = 3L
				),
				fragment_size = tfd_independent(
					tfd_bernoulli(logits = y2),
					reinterpreted_batch_ndims = 1L
				),
				position = tfd_independent(
					tfd_bernoulli(logits = y3),
					reinterpreted_batch_ndims = 1L
				)
			)
		}
	})
} # decoder_model_vae_20191216b


#' decoder_model_vae_vplot_position
#'
decoder_model_vae_vplot_position <- function(
	input_dim, 
	feature_dim, 
	filters0 = 64, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 2L), 
	name = NULL
){

	input_dim0 <- input_dim / prod(input_strides)
	feature_dim0 <- feature_dim / prod(feature_strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(feature_dim0, input_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(input_strides[1], feature_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			padding = 'same'
		)

		self$dense_2 <- layer_dense(units = feature_dim)
		self$bn_3 <- layer_batch_normalization()
		self$dense_3 <- layer_dense(units = feature_dim)

		function(x, mask = NULL){

			y1 <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3()

			y2 <- x %>%
				self$dense_2() %>%
				self$bn_3()

			list(
				vplot =	tfd_independent(
					tfd_bernoulli(logits = y1), 
					reinterpreted_batch_ndims = 3L
				),
				position = tfd_independent(
					tfd_bernoulli(logits = y2),
					reinterpreted_batch_ndims = 1L
				)
			)
		}
	})
} # decoder_model_vae_vplot_position


#' decoder_model_vae_position_fragment_size
#'
decoder_model_vae_position_fragment_size <- function(
	input_dim, 
	feature_dim, 
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$dense_2 <- layer_dense(units = 8L, activation = 'relu')
		self$dense_3 <- layer_dense(units = feature_dim)

		self$dense_4 <- layer_dense(units = 8L, activation = 'relu')
		self$dense_5 <- layer_dense(units = feature_dim)

		function(x, mask = NULL){

			y2 <- x %>%
				self$dense_2() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_3()

			y3 <- x %>%
				self$dense_4() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_5()

			list(
				fragment_size = tfd_independent(
					tfd_bernoulli(logits = y2),
					reinterpreted_batch_ndims = 1L
				),
				position = tfd_independent(
					tfd_bernoulli(logits = y3),
					reinterpreted_batch_ndims = 1L
				)
			)
		}
	})
} # decoder_model_vae_position_fragment_size


