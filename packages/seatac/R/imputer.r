#' imputer_model_baseline
#'
imputer_model_baseline <- function(
	input_dim,
	feature_dim,
	hidden_dim = 64L,
	name = NULL
){

	  keras_model_custom(name = name, function(self){

			self$dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')
			self$dense_2 <- layer_dense(units = input_dim * feature_dim, activation = 'relu')
			self$reshape_1 <- layer_reshape(target_shape = c(feature_dim, input_dim, 1L))

			function(x, mask = NULL){

				y <- x %>%
					layer_flatten() %>%
					self$dense_1() %>%
					layer_dropout(rate = 0.2) %>%
					self$dense_2() %>%
					layer_dropout(rate = 0.2) %>%
					self$reshape_1()
				y
			}
	})
} # decoder_model_vae_position_fragment_size

#' imputer_model_fragment_size_position
#'
imputer_model_fragment_size_position <- function(
	input_dim,
	feature_dim,
	hidden_dim = 8L,
	name = NULL
){

  keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$dense_2 <- layer_dense(units = feature_dim)
		self$dense_3 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$dense_4 <- layer_dense(units = input_dim)

		function(x, mask = NULL){

			y1 <- x[[1]] %>%
				self$dense_1() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_2() %>%
				layer_reshape(target_shape = c(feature_dim, 1L, 1L))

			y2 <- x[[2]] %>%
				self$dense_3() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_4() %>%
				layer_reshape(target_shape = c(1L, input_dim, 1L))

			y <- (y1 + y2) %>%
				layer_activation(activation = 'relu')

			y
		}
	})
} # imputer_model_fragment_size_position

#' imputer_model_position_fragment_size
#'
imputer_model_position_fragment_size <- function(
	input_dim,
	feature_dim,
	hidden_dim = 8L,
	name = NULL
){

  keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$dense_2 <- layer_dense(units = feature_dim)
		self$dense_3 <- layer_dense(units = hidden_dim, activation = 'relu')
		self$dense_4 <- layer_dense(units = input_dim)

		function(x, mask = NULL){

			browser()

			y1 <- x[[1]] %>%
				self$dense_1() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_2() %>%
				layer_reshape(target_shape = c(feature_dim, 1L, 1L))

			y2 <- x[[2]] %>%
				self$dense_3() %>%
				layer_dropout(rate = 0.2) %>%
				self$dense_4() %>%
				layer_reshape(target_shape = c(1L, input_dim, 1L))

			y <- (y1 + y2) %>%
				layer_activation(activation = 'relu')

			y
		}
	})
} # imputer_model_position_fragment_size

