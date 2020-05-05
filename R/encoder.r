#' encoder model for the vplot
#'
encoder_model <- function(
	latent_dim, 
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	window_strides = c(2L, 2L, 2L), 
	interval_strides = c(2L, 2L, 1L), 
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = c(window_strides[1], interval_strides[1]),
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$conv_2 <- layer_conv_2d(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(window_strides[2], interval_strides[2]),
			activation = 'relu'
		)

		self$bn_2 <- layer_batch_normalization()

		self$conv_3 <- layer_conv_2d(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(window_strides[3], interval_strides[3]),
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

