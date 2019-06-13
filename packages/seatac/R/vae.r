vae <- function(M, latent_dims = 10){

  sampling <- function(arg){
    z_mean <- arg[, 1:(latent_dims)]
    z_log_var <- arg[, (latent_dims+ 1):(2 * latent_dims)]
    epsilon <- k_random_normal(
      shape = c(k_shape(z_mean)[[1]]), 
      mean = 0.,
      stddev = 1.0
    )
    z_mean + k_exp(z_log_var / 2) * epsilon
  }

	input_layer <- layer_input(shape = shape(M))
	encoder_layers <- input_layer %>% 
		layer_reshape(target_shape = shape(M, 1)) %>% 
		layer_conv_1d(filters = 5, kernel_size = 10,  activation = 'relu', padding = 'same') %>%
		layer_max_pooling_1d(pool_size = 4) %>% 
		layer_flatten() 

	z_mean <- layer_dense(encoder_layers, latent_dims)
	z_log_var <- layer_dense(encoder_layers, latent_dims)
	z <- layer_concatenate(list(z_mean, z_log_var)) %>% layer_lambda(sampling)
		decoder_layers <- layer_dense(units = M, activation = 'relu')
		output_layer <- decoder_layers(z)

	vae <- keras_model(input_layer, output_layer)

	vae_loss <- function(x, x_decoded_mean){
		xent_loss <- (M / 1.0) * loss_mean_squared_error(x, x_decoded_mean) # the reconstruction loss
		kl_loss <- -0.5 * k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
		xent_loss + kl_loss
	}
	vae %>% compile(optimizer = 'rmsprop', loss = vae_loss)

} #  vae
