#' cvae_encoder_model
#'
cvae_encoder_model <- function(x, output_dim = 16L){

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

} # cvae_encoder_model


#' encoding model for the vplot
#'
gmm_cvae_encoder_model <- function(
	latent_dim, 
	sequence_dim,
	window_size,
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	name = NULL
){

		keras_model_custom(name = name, function(self){

			self$conv_1 <- layer_conv_2d(
				filters = filters[1],
				kernel_size = kernel_size[1],
				strides = shape(input_strides[1], feature_strides[1]),
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

			self$dense_1 <- layer_dense(units = sequence_dim)

			self$dense_2 <- layer_dense(units = 2 * latent_dim)

			self$sequence_conv <- layer_conv_1d(
				filters = 300L,
				kernel_size = 26L,
				strides = 1L,
				activation = 'relu'
			)

			self$bigru <- bidirectional(layer = layer_cudnn_gru(units = sequence_dim / 2))

			self$dense_3 <- layer_dense(units = latent_dim)

			function(x, mask = NULL){

				h_vplot <- x[[1]] %>% 
					self$conv_1() %>%
					self$bn_1() %>%
					self$conv_2() %>%
					self$bn_2() %>%
					self$conv_3() %>%
					self$bn_3() %>%
					layer_flatten() %>%
					self$dense_1()

				h_sequence <- x[[2]] %>% 
					layer_embedding(
						input_dim = 5L,
						output_dim = 5L,
						input_length = window_size,
						weights = list(diag(5L)),
						trainable = FALSE
					)  %>%
					self$sequence_conv() %>%
					layer_max_pooling_1d(pool_size = 13L, strides = 13L) %>%
					layer_dropout(0.2) %>%
					self$bigru() %>%
					layer_dropout(0.5)

				y <- layer_add(list(h_vplot, h_sequence)) %>%
					self$dense_2()

				z_sequence <- h_sequence %>%
					self$dense_3()

				list(
					z = tfd_multivariate_normal_diag(
						loc = y[, 1:latent_dim],
						scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
					),
					z_sequence = z_sequence
				)
			}
	})
}

