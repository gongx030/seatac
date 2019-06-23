#' build_model
#'

build_model <- function(input_dim, feature_dim, gpu = TRUE){
	model_20190621a(input_dim, feature_dim, gpu)
}

model_20190621a <- function(input_dim, feature_dim, gpu = TRUE){

	input_layer <- layer_input(shape = shape(input_dim, feature_dim))
  conv_1d_output_layer <- input_layer %>% 
  		layer_reshape(target_shape = c(input_dim, feature_dim, 1)) %>%
      time_distributed(layer_conv_1d(filters = 4L, kernel_size = 10L, strides = 1L, activation = 'relu')) %>%
      time_distributed(layer_flatten())

  if (gpu){
    flog.info('use GPU for training')
  	gru_output_layer <- conv_1d_output_layer %>%
#      bidirectional(layer_cudnn_gru(units = 16L, return_sequences = TRUE)) %>%
#      bidirectional(layer_cudnn_gru(units = 8L, return_sequences = TRUE)) %>%
      bidirectional(layer_cudnn_gru(units = 16L, return_sequences = TRUE))
  }else{
#    flog.info('use CPU for training')
#  	gru_output_layer <- conv_1d_output_layer %>%
#      bidirectional(layer_gru(units = gru_units, return_sequences = TRUE)) %>%
#      bidirectional(layer_gru(units = gru_units, return_sequences = TRUE)) %>%
#      bidirectional(layer_gru(units = gru_units, return_sequences = TRUE))
  }


#  left_layer <- gru_output_layer %>% 
#      time_distributed(layer_dense(units = feature_dim, activation = 'softmax')) 

  center_layer <- gru_output_layer %>% 
      time_distributed(layer_dense(units = feature_dim, activation = 'softmax')) 

#  right_layer <- gru_output_layer %>% 
#      time_distributed(layer_dense(units = feature_dim, activation = 'softmax')) 

#  model <- keras_model(input_layer, list(left_layer, center_layer, right_layer))
  model <- keras_model(input_layer, center_layer)
	print(model)
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
