#' VplotEncoder
#'
#' A Vplot encoder network
#' 
#' @param filters filters
#' @param kernel_size kernel size
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param distribution Output distributions (MultivariateNormalDiag, LogNormal, or None)
#' @param rate Dropout rate (default: 0.1)
#' @param name model name
#'
VplotEncoder <- function(
	downsample_layers = 4L,
	filters = 32L,
	kernel_size = 3L,
	rate = 0.1,
	momentum = 0.8,
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$downsample_layers <- downsample_layers

		self$conv <- lapply(1:self$downsample_layers, function(i) tf$keras$Sequential(list(
			tf$keras$layers$Conv2D(
				filters = filters,
				kernel_size = kernel_size,
				strides = shape(2L, 2L)
			),
			tf$keras$layers$Activation('relu'),
			tf$keras$layers$BatchNormalization(momentum = momentum)
		)))

		self$flatten_1 <- tf$keras$layers$Flatten()

		function(x, ...){

			for (i in 1:self$downsample_layers){
				x <- x %>% 
					self$conv[[i - 1]]() 
			}
			y <- x %>% 
				self$flatten_1()

			y
				
		}
	})
}


