#' composite autoencoder
#'
cae <- function(input_dim, feature_dim, latent_dim = 10){

	gru_units <- 16L

	input_layer <- layer_input(shape = shape(input_dim, feature_dim))
	output_layer <- input_layer %>% 
		layer_reshape(target_shape = c(input_dim, feature_dim, 1)) %>%
		time_distributed(layer_conv_1d(filters = 16, kernel_size = 5, activation = 'relu')) %>%
		time_distributed(layer_conv_1d(filters = 8, kernel_size = 2, activation = 'relu')) %>%
		time_distributed(layer_flatten()) %>%
    bidirectional(layer_cudnn_gru(units = gru_units, return_sequences = TRUE)) %>%
    bidirectional(layer_cudnn_gru(units = gru_units, return_sequences = TRUE)) %>%
    bidirectional(layer_cudnn_gru(units = gru_units, return_sequences = TRUE)) %>%
    time_distributed(layer_dense(units = feature_dim, activation = 'softmax'))

  model <- keras_model(input_layer, output_layer)
	print(model)
#	model %>% compile(optimizer = 'adam', loss = weighted_binary_crossentropy)
	model %>% compile(optimizer = 'adam', loss = 'binary_crossentropy')
  model

} #  vae

weighted_binary_crossentropy <- function(y_true, y_pred){
  one_weight <- 1
  zero_weight <- 0.1
  K <- backend()
  b_ce <- K$binary_crossentropy(y_true, y_pred)
  weight_vector <- y_true * one_weight  + (1. - y_true) * zero_weight
  weighted_b_ce <- weight_vector * b_ce
  K$mean(weighted_b_ce)
}
