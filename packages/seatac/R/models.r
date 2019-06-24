#' build_model
#'

build_model <- function(input_dim, feature_dim, gpu = TRUE){
	model_20190621a(input_dim, feature_dim, gpu)
}

model_20190621a <- function(input_dim, feature_dim, gpu = TRUE){

	input_layer <- layer_input(shape = shape(input_dim, feature_dim))
  conv_1d_output_layer <- input_layer %>% 
  		layer_reshape(target_shape = c(input_dim, feature_dim, 1)) %>%
      time_distributed(layer_conv_1d(filters = 4L, kernel_size = 5L, strides = 1L, activation = 'relu')) %>%
      time_distributed(layer_flatten()) %>%
      time_distributed(layer_dropout(rate = 0.3))

  if (gpu){
    flog.info('use GPU for training')
  	gru_output_layer <- conv_1d_output_layer %>%
      bidirectional(layer_cudnn_gru(units = 16L, return_sequences = TRUE))
  }else{
    flog.info('use CPU for training')
  	gru_output_layer <- conv_1d_output_layer %>%
      bidirectional(layer_gru(units = 16L, return_sequences = TRUE)) 
  }

  left_layer <- gru_output_layer %>% 
      time_distributed(layer_dropout(rate = 0.3)) %>%
      time_distributed(layer_dense(units = feature_dim, activation = 'softmax')) 

  center_layer <- gru_output_layer %>% 
      time_distributed(layer_dropout(rate = 0.3)) %>%
      time_distributed(layer_dense(units = feature_dim, activation = 'softmax')) 

  right_layer <- gru_output_layer %>% 
      time_distributed(layer_dropout(rate = 0.3)) %>%
      time_distributed(layer_dense(units = feature_dim, activation = 'softmax')) 

  model <- keras_model(input_layer, list(left_layer, center_layer, right_layer))
	print(model)
	model %>% compile(optimizer = 'adam', loss = c(loss_binary_crossentropy, loss_binary_crossentropy, loss_binary_crossentropy))
  model

} #  model_20190621a

weighted_binary_crossentropy <- function(y_true, y_pred){
  one_weight <- 1
  zero_weight <- 0.1
  K <- backend()
  b_ce <- K$binary_crossentropy(y_true, y_pred)
  weight_vector <- y_true * one_weight  + (1. - y_true) * zero_weight
  weighted_b_ce <- weight_vector * b_ce
  K$mean(weighted_b_ce)
}


vanilla_vae <- function(input_dim, feature_dim, latent_dim){

	encoder_model <- function(latent_dim, name = NULL){
		keras_model_custom(name = name, function(self){
			self$conv_1 <- layer_conv_2d(
				filters = 4L,
				kernel_size = 3L,
				strides = 2L,
				activation = 'relu'
			)
			self$flatten_1 <- layer_flatten()
			self$dense_1 <- layer_dense(units = 2 * latent_dim)
			function(x, mask = NULL){
				x <- x %>% 
					self$conv_1() %>%
					self$flatten_1() %>%
					self$dense_1()
				tfd_multivariate_normal_diag(
					loc = x[, 1:latent_dim],
					scale_diag = tf$nn$softplus(x[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
				)
			}
		})
  }
  encoder <- encoder_model(latent_dim)

	decoder_model <- function(input_dim, feature_dim, name = NULL){

		keras_model_custom(name = name, function(self){
			self$dense_1 <- layer_dense(units = 8 * 8 * 16, activation = "relu")
			self$reshape_1 <- layer_reshape(target_shape = c(8L, 8L, 16L))
			self$deconv_1 <- layer_conv_2d_transpose(
				filters = 16L,
				kernel_size = 3L,
				strides = 4L,
				padding = 'same',
				activation = 'relu'
			)
			self$deconv_3 <- layer_conv_2d_transpose(
				filters = 1L,
				kernel_size = 3L,
				strides = 1L,
				padding = 'same'
			)
			function(x, mask = NULL){
				x <- x %>% 
					self$dense_1() %>%
					self$reshape_1() %>%
					self$deconv_1() %>%
					self$deconv_3()
				
				tfd_independent(
					tfd_bernoulli(logits = x), 
					reinterpreted_batch_ndims = 3L
				)
			}
		})
	}
	decoder <- decoder_model(input_dim, feature_dim)
	latent_prior <- tfd_multivariate_normal_diag(loc = rep(0, latent_dim), scale_identity_multiplier = 1)

	list(encoder = encoder, decoder = decoder, latent_prior = latent_prior, rainable_prior = FALSE)

} # vanilla_vae


gmm_vae <- function(input_dim, feature_dim, latent_dim, n_components){

	learnable_prior_model <- function(latent_dim, n_components, name = NULL){
		keras_model_custom(name = name, function(self) {
			self$loc <- tf$get_variable(
				name = 'loc',
				shape = list(n_components, latent_dim),
				dtype = tf$float32
			)
			self$raw_scale_diag <- tf$get_variable(
				name = 'raw_scale_diag',
				shape = list(n_components, latent_dim),
				dtype = tf$float32
			)
			self$mixture_logits <- tf$get_variable(
				name = 'mixture_logits',
				shape = c(n_components),
				dtype = tf$float32
			)
			function (x, mask = NULL) {
				tfd_mixture_same_family(
					components_distribution = tfd_multivariate_normal_diag(loc = self$loc, scale_diag = tf$nn$softplus(self$raw_scale_diag)),
					mixture_distribution = tfd_categorical(logits = self$mixture_logits)
				)
			}
		})
	}
	latent_prior_model <- learnable_prior_model(latent_dim, n_components)

	encoder_model <- function(latent_dim, name = NULL){
		keras_model_custom(name = name, function(self){
			self$conv_1 <- layer_conv_2d(
				filters = 4L,
				kernel_size = 3L,
				strides = 2L,
				activation = 'relu'
			)
			self$flatten_1 <- layer_flatten()
			self$dense_1 <- layer_dense(units = 2 * latent_dim)
			function(x, mask = NULL){
				x <- x %>% 
					self$conv_1() %>%
					self$flatten_1() %>%
					self$dense_1()
				tfd_multivariate_normal_diag(
					loc = x[, 1:latent_dim],
					scale_diag = tf$nn$softplus(x[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
				)
			}
		})
  }
  encoder <- encoder_model(latent_dim)

	decoder_model <- function(input_dim, feature_dim, name = NULL){

		keras_model_custom(name = name, function(self){
			self$dense_1 <- layer_dense(units = 8 * 8 * 16, activation = "relu")
			self$reshape_1 <- layer_reshape(target_shape = c(8L, 8L, 16L))
			self$deconv_1 <- layer_conv_2d_transpose(
				filters = 16L,
				kernel_size = 3L,
				strides = 4L,
				padding = 'same',
				activation = 'relu'
			)
			self$deconv_3 <- layer_conv_2d_transpose(
				filters = 1L,
				kernel_size = 3L,
				strides = 1L,
				padding = 'same'
			)
			function(x, mask = NULL){
				x <- x %>% 
					self$dense_1() %>%
					self$reshape_1() %>%
					self$deconv_1() %>%
					self$deconv_3()
				
				tfd_independent(
					tfd_bernoulli(logits = x), 
					reinterpreted_batch_ndims = 3L
				)
			}
		})
	}
	decoder <- decoder_model(input_dim, feature_dim)

	list(encoder = encoder, decoder = decoder, latent_prior_model = latent_prior_model, trainable_prior = TRUE)

} # gmm_vae

