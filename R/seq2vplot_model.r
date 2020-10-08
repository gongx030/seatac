#' Seq2VplotModel
#'

Seq2VplotModel <- function(
	d_model,
	num_layers, 
	num_heads, 
	dff, 
	kmers_size,
	vae,
	step_size = 2L,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$d_model <- d_model
		stopifnot(self$d_model %% num_heads == 0)
		self$num_layers <- num_layers
		self$vae <- vae
		self$vae$trainable <- FALSE
		self$pos_encoding <- positional_encoding(self$vae$block_size, self$d_model) %>%
		tf$cast(tf$float32) %>%
		tf$expand_dims(0L)
		self$embedding <- tf$keras$layers$Embedding(as.integer(kmers_size), self$d_model)
		self$dropout_1 <- tf$keras$layers$Dropout(rate)
	
		if (num_layers > 0){
			self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(self$d_model, num_heads = num_heads, dff = dff, rate = rate))
		}
	
		self$dense_1 <- tf$keras$layers$Dense(units = self$vae$decoder$vplot_decoder$vplot_height, activation = 'relu')
		self$pool <- tf$keras$layers$MaxPooling1D(self$vae$bin_size, self$vae$bin_size)
	
		function(x, training = TRUE){
			x <- x %>% self$embedding()
			x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))
			x <- x + self$pos_encoding
			x <- x %>% self$dropout_1()
	
			if (num_layers > 0){
				for (i in seq_len(self$num_layers)){
					x <- self$enc_layers[[i - 1]](x)	# zero-based
				}
			}

			x <- x %>% self$dense_1()
			x <- x %>% self$pool()
			x <- x %>% tf$transpose(shape(0L, 2L, 1L))
			x <- x %>% tf$expand_dims(3L)
			posterior <- self$vae$encoder(x)
			z <- posterior$mean()
			res <- self$vae$decoder(z)
			res$z <- z
			res
		}	
	})
}
