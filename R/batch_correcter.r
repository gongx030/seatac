#' batch correcter
#'
#' batch correcter for the V-plot
#' 
#' Adopted from the batch correction layer used in scScope
#'
batch_correcter_model <- function(
	name = NULL,
	n_samples = NULL,
	window_dim = NULL,
	interval_dim = NULL
){
	keras_model_custom(name = name, function(self){

		self$loc <- tf$Variable(tf$random$normal(shape(n_samples, 1L, interval_dim)))

		function(x, mask = NULL, training = TRUE){

			y <- self$loc %>%  
				tf$tile(shape(1L, window_dim, 1L)) %>%
				tf$reshape(shape(n_samples, window_dim * interval_dim))

			y <- tf$matmul(x[[2]], y) %>%
				layer_reshape(shape(window_dim, interval_dim, 1L))

			y <- (x[[1]] - y) %>%
				layer_activation(activation = 'relu')
			
			y

		}
	})

} # batch_correcter_model
