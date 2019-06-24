vae <- function(input_dim, feature_dim, latent_dim){

	encoder_model <- function(latent_dim, name = NULL){
		keras_model_custom(name = name, function(self){
			self$conv_1 <- layer_conv_2d(
				filters = 32,
				kernel_size = c(3, 3),
				activation = 'relu'
			)
			self$dense <- layer_dense(units = 2 * latent_dim)
			function(x, mask = NULL){
				x <- x %>% 
					self$conv_1() %>%
					layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
					layer_dropout(rate = 0.25) %>%
					layer_flatten() %>%
					self$dense()
				tfd_multivariate_normal_diag(
					loc = x[, 1:latent_dim],
					scale_diag = tf$nn$softplus(x[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
				)
			}
		})
  }
  encoder <- encoder_model(latent_dim)
	browser()
} # vae
