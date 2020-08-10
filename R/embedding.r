kmer_embedding_model <- function(
	kmer_dim,
	bin_size,
	input_dim,
	context_dim,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$embedding <- layer_embedding(
			input_dim = input_dim,
			output_dim = kmer_dim
		)

#		self$bigru <- bidirectional(layer = layer_cudnn_gru(
#			units = context_dim,
#			return_sequence = TRUE
#		))

		function(x, mask = NULL, training = TRUE){

			y <- x %>%
				self$embedding() %>%
				layer_average_pooling_1d(pool_size = bin_size, strides = bin_size)
#				self$bigru()
			y

		}
	})

} # encoder_model


