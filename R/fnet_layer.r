#' FNetLayer
#'
#' A FNet encoder layer
#' 
#' @param embedding_dim Embedding dimension (default: 256L)
#' @param rate Dropout rate (default: 0.1)
#' @param name model name
#'
FNetLayer <- function(
	embedding_dim = 256L,	
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self) {
		
		self$ffn <- tf$keras$Sequential(list(
			tf$keras$layers$Dense(units = embedding_dim),
			tf$keras$layers$Activation('gelu'),
			tf$keras$layers$Dropout(rate = rate),
      tf$keras$layers$Dense(units = embedding_dim)
		))

		self$normalize1 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
    self$normalize2 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)

		function(x, ...){

			xd <- x %>% 
				tf$cast(dtype = tf$dtypes$complex64) %>%
				tf$signal$fft2d() %>%
				tf$cast(dtype = tf$dtypes$float32)
			
			x <- x + xd

			x <- x %>% self$normalize1()

			x_ffn <- x %>% self$ffn()

			x <- x + x_ffn

			x <- x %>% self$normalize2()

			x

		}
	})
}


