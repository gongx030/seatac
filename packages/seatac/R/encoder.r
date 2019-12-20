#' encoder model for the vplot
#'
encoder_model_cvae <- function(
	latent_dim, 
	window_size,
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	hidden_dim = 8L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = c(input_strides[1], feature_strides[1]),
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$conv_2 <- layer_conv_2d(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			activation = 'relu'
		)

		self$bn_2 <- layer_batch_normalization()

		self$conv_3 <- layer_conv_2d(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			activation = 'relu'
		)

		self$bn_3 <- layer_batch_normalization()

		self$dense_1 <- layer_dense(units = hidden_dim)

		self$dense_2 <- layer_dense(units = hidden_dim, activation = 'softmax')

		self$dense_3 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			h_vplot  <- x[[1]] %>% 
				self$conv_1() %>%
				self$bn_1(training = training) %>%
				self$conv_2() %>%
				self$bn_2(training = training) %>%
				self$conv_3() %>%
				self$bn_3(training = training) %>%
				layer_flatten() %>%
				self$dense_1()

			h_fragment_size <- x[[2]] %>%
				self$dense_2()

			y <- layer_concatenate(list(h_vplot, h_fragment_size)) %>%
				self$dense_3()

			list(
				z = tfd_multivariate_normal_diag(
					loc = y[, 1:latent_dim],
					scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
				),
				fragment_size = h_fragment_size
			)
		}
	})

} # encoder_model_cvae


#' encoder model for the vplot
#'
encoder_model_vae <- function(
	latent_dim, 
	window_size,
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	hidden_dim = 8L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = c(input_strides[1], feature_strides[1]),
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$conv_2 <- layer_conv_2d(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			activation = 'relu'
		)

		self$bn_2 <- layer_batch_normalization()

		self$conv_3 <- layer_conv_2d(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			activation = 'relu'
		)

		self$bn_3 <- layer_batch_normalization()

		self$dense_1 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y <- x %>% 
				self$conv_1() %>%
				self$bn_1(training = training) %>%
				self$conv_2() %>%
				self$bn_2(training = training) %>%
				self$conv_3() %>%
				self$bn_3(training = training) %>%
				layer_flatten() %>%
				self$dense_1()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})

} # encoder_model_vae


#' encoder model for the vplot
#'
encoder_model_vae_baseline <- function(
	latent_dim, 
	window_size,
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	hidden_dim = 8L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = c(input_strides[1], feature_strides[1]),
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$conv_2 <- layer_conv_2d(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			activation = 'relu'
		)

		self$bn_2 <- layer_batch_normalization()

		self$conv_3 <- layer_conv_2d(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			activation = 'relu'
		)

		self$bn_3 <- layer_batch_normalization()

		self$dense_1 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y <- x %>% 
				self$conv_1() %>%
				self$bn_1(training = training) %>%
				self$conv_2() %>%
				self$bn_2(training = training) %>%
				self$conv_3() %>%
				self$bn_3(training = training) %>%
				layer_flatten() %>%
				self$dense_1()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-3)
			)
		}
	})

} # encoder_model_vae_baseline


#' encoder model for the vplot
#'
encoder_model_vae_20191216a <- function(
	latent_dim, 
	window_size,
	hidden_dim = 8L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$bn_1 <- layer_batch_normalization()

		self$dense_2 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$bn_2 <- layer_batch_normalization()

		self$dense_3 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y1 <- x[[1]] %>% 
				self$dense_1() %>%
				self$bn_1()

			y2 <- x[[2]] %>% 
				self$dense_2() %>%
				self$bn_2()

			y <- layer_concatenate(list(y1, y2)) %>%
				self$dense_3()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})

} # encoder_model_vae_20191216a


#' encoder model for the vplot
#'
encoder_model_vae_20191216b <- function(
	latent_dim, 
	window_size,
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	hidden_dim = 8L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = c(input_strides[1], feature_strides[1]),
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$conv_2 <- layer_conv_2d(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			activation = 'relu'
		)

		self$bn_2 <- layer_batch_normalization()

		self$conv_3 <- layer_conv_2d(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			activation = 'relu'
		)

		self$bn_3 <- layer_batch_normalization()

		self$dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$bn_4 <- layer_batch_normalization()

		self$dense_2 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$bn_5 <- layer_batch_normalization()

		self$dense_3 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y1 <- x[[1]] %>% 
				self$conv_1() %>%
				self$bn_1(training = training) %>%
				self$conv_2() %>%
				self$bn_2(training = training) %>%
				self$conv_3() %>%
				self$bn_3(training = training) %>%
				layer_flatten() 
			
			y2 <- x[[2]] %>%
				self$dense_1() %>%
				self$bn_4()

			y3 <- x[[3]] %>%
				self$dense_2() %>%
				self$bn_5()

			y <- layer_concatenate(list(y1, y2, y3)) %>%
				self$dense_3()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}	
	})

} # encoder_model_vae_20191216b


#' encoder_model_vae_position_fragment_size
#'
encoder_model_vae_position_fragment_size <- function(
	latent_dim, 
	window_size,
	hidden_dim = 32L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$bn_1 <- layer_batch_normalization()

		self$dense_2 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$bn_2 <- layer_batch_normalization()

		self$dense_3 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y1 <- x[[1]] %>% 
				self$dense_1() 

			y2 <- x[[2]] %>% 
				self$dense_2() 

			y <- layer_concatenate(list(y1, y2)) %>%
				self$dense_3()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})

} # encoder_model_vae_position_fragment_size


#' encoder_model_vae_position_fragment_size_cnn
#'
encoder_model_vae_position_fragment_size_cnn <- function(
	latent_dim, 
	fragment_size_filters = c(32L, 32L, 32L),
	fragment_size_kernel_sizes = c(3L, 3L, 3L),
	fragment_size_strides = c(2L, 2L, 2L),
	position_filters = c(32L, 32L, 32L),
	position_kernel_sizes = c(3L, 3L, 3L),
	position_strides = c(2L, 2L, 2L),
	hidden_dim = 32L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$fragment_size_conv_1 <- layer_conv_1d(
			filters = fragment_size_filters[1],
			kernel_size = fragment_size_kernel_sizes[1],
			strides = fragment_size_kernel_sizes[1],
			activation = 'relu'
		)

		self$fragment_size_bn_1 <- layer_batch_normalization()

		self$fragment_size_conv_2 <- layer_conv_1d(
			filters = fragment_size_filters[2],
			kernel_size = fragment_size_kernel_sizes[2],
			strides = fragment_size_strides[2],
			activation = 'relu'
		)

		self$fragment_size_bn_2 <- layer_batch_normalization()

		self$fragment_size_conv_3 <- layer_conv_1d(
			filters = fragment_size_filters[3],
			kernel_size = fragment_size_kernel_sizes[3],
			strides = fragment_size_strides[3],
			activation = 'relu'
		)

		self$fragment_size_bn_3 <- layer_batch_normalization()

		self$fragment_size_dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')

		self$position_conv_1 <- layer_conv_1d(
			filters = position_filters[1],
			kernel_size = position_kernel_sizes[1],
			strides = position_strides[1],
			activation = 'relu'
		)

		self$position_bn_1 <- layer_batch_normalization()

		self$position_conv_2 <- layer_conv_1d(
			filters = position_filters[2],
			kernel_size = position_kernel_sizes[2],
			strides = position_strides[2],
			activation = 'relu'
		)

		self$position_bn_2 <- layer_batch_normalization()

		self$position_conv_3 <- layer_conv_1d(
			filters = position_filters[3],
			kernel_size = position_kernel_sizes[3],
			strides = position_strides[3],
			activation = 'relu'
		)

		self$position_bn_3 <- layer_batch_normalization()

		self$position_dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')

		self$merge_dense_1 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y1 <- x[[1]] %>% 
				k_expand_dims() %>% 
				self$fragment_size_conv_1() %>%
				self$fragment_size_bn_1(training = training) %>%
				self$fragment_size_conv_2() %>%
				self$fragment_size_bn_2(training = training) %>%
				self$fragment_size_conv_3() %>%
				self$fragment_size_bn_3(training = training) %>%
				layer_flatten() %>%
				self$fragment_size_dense_1()

			y2 <- x[[2]] %>% 
				k_expand_dims() %>% 
				self$position_conv_1() %>%
				self$position_bn_1(training = training) %>%
				self$position_conv_2() %>%
				self$position_bn_2(training = training) %>%
				self$position_conv_3() %>%
				self$position_bn_3(training = training) %>%
				layer_flatten() %>%
				self$position_dense_1()

			y <- list(y1, y2) %>% 
				layer_concatenate() %>%
				layer_dropout(rate = 0.2) %>%
				self$merge_dense_1()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})

} # encoder_model_vae_position_fragment_size_cnn



#' encoder_model_vae_fragment_size_cnn
#'
encoder_model_vae_fragment_size_cnn <- function(
	latent_dim, 
	fragment_size_filters = c(32L, 32L, 32L),
	fragment_size_kernel_sizes = c(3L, 3L, 3L),
	fragment_size_strides = c(2L, 2L, 2L),
	hidden_dim = 32L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$fragment_size_conv_1 <- layer_conv_1d(
			filters = fragment_size_filters[1],
			kernel_size = fragment_size_kernel_sizes[1],
			strides = fragment_size_kernel_sizes[1],
			activation = 'relu'
		)

		self$fragment_size_bn_1 <- layer_batch_normalization()

		self$fragment_size_conv_2 <- layer_conv_1d(
			filters = fragment_size_filters[2],
			kernel_size = fragment_size_kernel_sizes[2],
			strides = fragment_size_strides[2],
			activation = 'relu'
		)

		self$fragment_size_bn_2 <- layer_batch_normalization()

		self$fragment_size_conv_3 <- layer_conv_1d(
			filters = fragment_size_filters[3],
			kernel_size = fragment_size_kernel_sizes[3],
			strides = fragment_size_strides[3],
			activation = 'relu'
		)

		self$fragment_size_bn_3 <- layer_batch_normalization()

		self$fragment_size_dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')

		self$merge_dense_1 <- layer_dense(units = latent_dim * 2)

		function(x, mask = NULL, training = TRUE){

			y <- x %>% 
				k_expand_dims() %>% 
				self$fragment_size_conv_1() %>%
				self$fragment_size_bn_1(training = training) %>%
				self$fragment_size_conv_2() %>%
				self$fragment_size_bn_2(training = training) %>%
				self$fragment_size_conv_3() %>%
				self$fragment_size_bn_3(training = training) %>%
				layer_flatten() %>%
				self$fragment_size_dense_1() %>%
				layer_dropout(rate = 0.2) %>%
				self$merge_dense_1()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})

} # encoder_model_vae_fragment_size_cnn


