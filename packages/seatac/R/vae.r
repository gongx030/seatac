#' composite autoencoder
#'
cae <- function(input_dim, feature_dim, latent_dim = 10){

	input_layer <- layer_input(shape = shape(input_dim, feature_dim))
	encoder <- input_layer %>% 
    bidirectional(layer_gru(units = latent_dim, return_sequences = TRUE, activation = 'relu', dropout = 0.3)) %>%
    bidirectional(layer_gru(units = latent_dim, return_sequences = TRUE, activation = 'relu', dropout = 0.3))

  decoder <- encoder %>% 
    bidirectional(layer_gru(units = latent_dim, return_sequences = TRUE, activation = 'relu', dropout = 0.3)) %>%
    bidirectional(layer_gru(units = latent_dim, return_sequences = TRUE, activation = 'relu', dropout = 0.3)) %>%
#    time_distributed(layer_dense(units = 32, activation = 'relu')) %>%
    time_distributed(layer_dense(units = feature_dim, activation = 'softmax'))

  model <- keras_model(input_layer, decoder)
	model %>% compile(optimizer = 'adam', loss = weighted_binary_crossentropy)
#	model %>% compile(optimizer = 'adam', loss = 'binary_crossentropy')
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
