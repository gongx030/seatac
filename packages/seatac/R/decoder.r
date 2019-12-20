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


#' decoder_model_vae_baseline_conv
#'
decoder_model_vae_baseline_conv <- function(
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
} # decoder_model_vae_baseline_conv


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

		self$fragment_size_dense_1 <- layer_dense(units = 8L, activation = 'relu')
		self$fragment_size_dense_2 <- layer_dense(units = 32L, activation = 'relu')
		self$fragment_size_dense_3 <- layer_dense(units = feature_dim)

		self$position_dense_1 <- layer_dense(units = 8L, activation = 'relu')
		self$position_dense_2 <- layer_dense(units = 32L, activation = 'relu')
		self$position_dense_3 <- layer_dense(units = input_dim)

		function(x, mask = NULL){

			y1 <- x %>%
				self$fragment_size_dense_1() %>%
				layer_dropout(rate = 0.2) %>%
				self$fragment_size_dense_2() %>%
				layer_dropout(rate = 0.2) %>%
				self$fragment_size_dense_3()

			y2 <- x %>%
				self$position_dense_1() %>%
				layer_dropout(rate = 0.2) %>%
				self$position_dense_2() %>%
				layer_dropout(rate = 0.2) %>%
				self$position_dense_3()

			list(
				fragment_size = tfd_independent(
					tfd_bernoulli(logits = y1),
					reinterpreted_batch_ndims = 1L
				),
				position = tfd_independent(
					tfd_bernoulli(logits = y2),
					reinterpreted_batch_ndims = 1L
				)
			)
		}
	})
} # decoder_model_vae_position_fragment_size


#' decoder_model_vae_baseline_conv2
#'
#' The output layer is different from decoder_model_vae_baseline
#' This is used for weighting different entries in the loss function
#'
decoder_model_vae_baseline_conv2 <- function(
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

			tfd_bernoulli(logits = y)
		}
	})
} # decoder_model_vae_baseline_conv2


#' decoder_model_vae_baseline_dense
#'
decoder_model_vae_baseline_dense <- function(
	input_dim, 
	feature_dim, 
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = 32L, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$bn_1 <- layer_batch_normalization()
		
		self$dense_2 <- layer_dense(units = 128L, activation = 'relu')
		self$dropout_2 <- layer_dropout(rate = 0.2)
		self$bn_2 <- layer_batch_normalization()

		self$dense_3 <- layer_dense(units = input_dim * feature_dim, activation = 'relu')
		self$dropout_3 <- layer_dropout(rate = 0.2)
		self$bn_3 <- layer_batch_normalization()

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$bn_1() %>%

				self$dense_2() %>%
				self$dropout_2() %>%
				self$bn_2() %>%

				self$dense_3() %>%
				self$dropout_3() %>%
				self$bn_3() %>%

				layer_reshape(target_shape = c(feature_dim, input_dim, 1L))

			tfd_independent(
				tfd_bernoulli(logits = y), 
				reinterpreted_batch_ndims = 3L
			)
		}
	})
} # decoder_model_vae_baseline_dense


#' decoder_model_vae_vplot_parametric
#'
decoder_model_vae_vplot_parametric <- function(
	input_dim, 
	feature_dim, 
	n_components = 2L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$lambda <- tf$Variable(tf$random$normal(shape(1L)), name = 'lambda')

		self$dense_1 <- layer_dense(units = n_components, activation = 'softmax', name = 'mixture')

		self$fragment_size_index <- tf$constant(seq(0, 1, length.out = feature_dim))

		self$mono_mean <- tf$Variable(tf$random$normal(shape(1L)), name = 'mono_mean')
		self$mono_sd <- tf$Variable(tf$random$normal(shape(1L)), name = 'mono_sd')

		function(x, mask = NULL){

			probs <- x %>% self$dense_1()

			mono <- tf$exp(-tf$square(self$fragment_size_index - self$mono_mean) / tf$exp(self$mono_sd))

			nfr <- tf$exp(-tf$exp(self$lambda) * self$fragment_size_index)

			y <- tf$matmul(probs,  tf$stack(list(mono, nfr), axis = 1L), transpose_b = TRUE)

			tfd_independent(
				tfd_bernoulli(probs = y),
				reinterpreted_batch_ndims = 1L
			)
		}
	})
} # decoder_model_vae_vplot_parametric


#' decoder_model_vae_position_fragment_size_cnn
#'
decoder_model_vae_position_fragment_size_cnn <- function(

	input_dim, 
	feature_dim, 

	filters0 = 32L, 

	fragment_size_filters = c(32L, 32L, 1L), 
	fragment_size_kernel_sizes = c(3L, 3L, 3L), 
	fragment_size_strides = c(2L, 2L, 2L), 

	position_filters = c(32L, 32L, 1L), 
	position_kernel_sizes = c(3L, 3L, 3L), 
	position_strides = c(2L, 2L, 2L), 

	reinterpreted_batch_ndims = 1L,

	name = NULL
){

	fragment_size_dim0 <- feature_dim / prod(fragment_size_strides)
	position_dim0 <- input_dim / prod(position_strides)

	fragment_size_output_dim0 <- fragment_size_dim0 * filters0
	position_output_dim0 <- position_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$fragment_size_dense_1 <- layer_dense(units = fragment_size_output_dim0, activation = 'relu')

		self$fragment_size_deconv_1 <- layer_conv_2d_transpose(
			filters = fragment_size_filters[1],
			kernel_size = fragment_size_kernel_sizes[1],
			strides = shape(fragment_size_strides[1], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$fragment_size_bn_1 <- layer_batch_normalization()

		self$fragment_size_deconv_2 <- layer_conv_2d_transpose(
			filters = fragment_size_filters[2],
			kernel_size = fragment_size_kernel_sizes[2],
			strides = shape(fragment_size_strides[2], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$fragment_size_bn_2 <- layer_batch_normalization()

		self$fragment_size_deconv_3 <- layer_conv_2d_transpose(
			filters = fragment_size_filters[3],
			kernel_size = fragment_size_kernel_sizes[3],
			strides = shape(fragment_size_strides[3], 1L),
			padding = 'same'
		)

		self$position_dense_1 <- layer_dense(units = position_output_dim0, activation = 'relu')

		self$position_deconv_1 <- layer_conv_2d_transpose(
			filters = position_filters[1],
			kernel_size = position_kernel_sizes[1],
			strides = shape(position_strides[1], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$position_bn_1 <- layer_batch_normalization()

		self$position_deconv_2 <- layer_conv_2d_transpose(
			filters = position_filters[2],
			kernel_size = position_kernel_sizes[2],
			strides = shape(position_strides[2], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$position_bn_2 <- layer_batch_normalization()

		self$position_deconv_3 <- layer_conv_2d_transpose(
			filters = position_filters[3],
			kernel_size = position_kernel_sizes[3],
			strides = shape(position_strides[3], 1L),
			padding = 'same'
		)

		function(x, mask = NULL){

			y1 <- x %>%
				self$fragment_size_dense_1() %>%
				layer_dropout(0.2) %>%
				layer_reshape(target_shape = c(fragment_size_dim0, 1L, filters0)) %>%
				self$fragment_size_deconv_1() %>%
				self$fragment_size_bn_1() %>%
				self$fragment_size_deconv_2() %>%
				self$fragment_size_bn_2() %>%
				self$fragment_size_deconv_3() %>%
				layer_reshape(target_shape = c(feature_dim))

			y2 <- x %>%
				self$position_dense_1() %>%
				layer_dropout(0.2) %>%
				layer_reshape(target_shape = c(position_dim0, 1L, filters0)) %>%
				self$position_deconv_1() %>%
				self$position_bn_1() %>%
				self$position_deconv_2() %>%
				self$position_bn_2() %>%
				self$position_deconv_3() %>%
				layer_reshape(target_shape = c(input_dim))

			list(
				fragment_size = tfd_independent(
					tfd_bernoulli(logits = y1),
					reinterpreted_batch_ndims = reinterpreted_batch_ndims,
				),
				position = tfd_independent(
					tfd_bernoulli(logits = y2),
					reinterpreted_batch_ndims = reinterpreted_batch_ndims,
				)
			)
		}
	})
} # decoder_model_vae_position_fragment_size_cnn


#' decoder_model_vae_position_fragment_size_mixture_cnn
#'
decoder_model_vae_position_fragment_size_mixture_cnn <- function(

	input_dim, 
	feature_dim, 

	filters0 = 32L, 

	fragment_size_filters = c(32L, 32L, 1L), 
	fragment_size_kernel_sizes = c(3L, 3L, 3L), 
	fragment_size_strides = c(2L, 2L, 2L), 

	position_filters = c(32L, 32L, 1L), 
	position_kernel_sizes = c(3L, 3L, 3L), 
	position_strides = c(2L, 2L, 2L), 

	name = NULL
){

	position_dim0 <- input_dim / prod(position_strides)
	position_output_dim0 <- position_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$fragment_size_dense_1 <- layer_dense(units = 8L, activation = 'relu')
		self$fragment_size_dense_2 <- layer_dense(units = 32L, activation = 'relu')
		self$fragment_size_dense_3 <- layer_dense(units = feature_dim)

		self$position_dense_1 <- layer_dense(units = position_output_dim0, activation = 'relu')

		self$position_deconv_1 <- layer_conv_2d_transpose(
			filters = position_filters[1],
			kernel_size = position_kernel_sizes[1],
			strides = shape(position_strides[1], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$position_bn_1 <- layer_batch_normalization()

		self$position_deconv_2 <- layer_conv_2d_transpose(
			filters = position_filters[2],
			kernel_size = position_kernel_sizes[2],
			strides = shape(position_strides[2], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$position_bn_2 <- layer_batch_normalization()

		self$position_deconv_3 <- layer_conv_2d_transpose(
			filters = position_filters[3],
			kernel_size = position_kernel_sizes[3],
			strides = shape(position_strides[3], 1L),
			padding = 'same'
		)

		function(x, mask = NULL){

			y1 <- x %>%
				self$fragment_size_dense_1() %>%
				layer_dropout(rate = 0.2) %>%
				self$fragment_size_dense_2() %>%
				layer_dropout(rate = 0.2) %>%
				self$fragment_size_dense_3()

			y2 <- x %>%
				self$position_dense_1() %>%
				layer_dropout(0.2) %>%
				layer_reshape(target_shape = c(position_dim0, 1L, filters0)) %>%
				self$position_deconv_1() %>%
				self$position_bn_1() %>%
				self$position_deconv_2() %>%
				self$position_bn_2() %>%
				self$position_deconv_3() %>%
				layer_reshape(target_shape = c(input_dim))

			list(
				fragment_size = tfd_independent(
					tfd_bernoulli(logits = y1),
					reinterpreted_batch_ndims = 1L
				),
				position = tfd_independent(
					tfd_bernoulli(logits = y2),
					reinterpreted_batch_ndims = 1L
				)
			)
		}
	})
} # decoder_model_vae_position_fragment_size_mixture_cnn


#' decoder_model_vae_fragment_size_cnn
#'
decoder_model_vae_fragment_size_cnn <- function(

	feature_dim, 

	filters0 = 32L, 

	fragment_size_filters = c(32L, 32L, 1L), 
	fragment_size_kernel_sizes = c(3L, 3L, 3L), 
	fragment_size_strides = c(2L, 2L, 2L),

	name = NULL
){

	fragment_size_dim0 <- feature_dim / prod(fragment_size_strides)

	fragment_size_output_dim0 <- fragment_size_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$fragment_size_dense_1 <- layer_dense(units = fragment_size_output_dim0, activation = 'relu')

		self$fragment_size_deconv_1 <- layer_conv_2d_transpose(
			filters = fragment_size_filters[1],
			kernel_size = fragment_size_kernel_sizes[1],
			strides = shape(fragment_size_strides[1], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$fragment_size_bn_1 <- layer_batch_normalization()

		self$fragment_size_deconv_2 <- layer_conv_2d_transpose(
			filters = fragment_size_filters[2],
			kernel_size = fragment_size_kernel_sizes[2],
			strides = shape(fragment_size_strides[2], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$fragment_size_bn_2 <- layer_batch_normalization()

		self$fragment_size_deconv_3 <- layer_conv_2d_transpose(
			filters = fragment_size_filters[3],
			kernel_size = fragment_size_kernel_sizes[3],
			strides = shape(fragment_size_strides[3], 1L),
			padding = 'same'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$fragment_size_dense_1() %>%
				layer_dropout(0.2) %>%
				layer_reshape(target_shape = c(fragment_size_dim0, 1L, filters0)) %>%
				self$fragment_size_deconv_1() %>%
				self$fragment_size_bn_1() %>%
				self$fragment_size_deconv_2() %>%
				self$fragment_size_bn_2() %>%
				self$fragment_size_deconv_3() %>%
				layer_reshape(target_shape = c(feature_dim))

			fragment_size = tfd_independent(
				tfd_bernoulli(logits = y),
				reinterpreted_batch_ndims = 1L
			)
		}
	})
} # decoder_model_vae_fragment_size_cnn
