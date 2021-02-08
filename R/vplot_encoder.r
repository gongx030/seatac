#' VplotEncoder
#'
#' A Vplot encoder network
#' 
#' @param filters filters
#' @param kernel_size kernel size
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param name model name
#'
VplotEncoder <- function(
	filters = c(32L, 32L, 32L),
	kernel_size = c(3L, 3L, 3L),
	window_strides = c(2L, 2L, 2L),
	interval_strides = c(2L, 2L, 2L),
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$filters <- filters
		self$kernel_size <- kernel_size
		self$window_strides <-  window_strides
		self$interval_strides <-  interval_strides
		self$n_layers <- length(filters)

		self$conv1d <- lapply(1:self$n_layers, function(i) 
			tf$keras$layers$Conv2D(
				filters = filters[i],
				kernel_size = kernel_size[i],
				strides = shape(interval_strides[i], window_strides[i]),
				activation = 'relu'
			)
		)

		self$bn <- lapply(1:self$n_layers, function(i) tf$keras$layers$BatchNormalization())

		function(x, ...){

			for (i in 1:self$n_layers){
				x <- x %>% 
					self$conv1d[[i - 1]]() %>% # zero-based
					self$bn[[i - 1]]()	# zero-based
			}
			x
		}
	})
}


