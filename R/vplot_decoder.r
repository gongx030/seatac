#' VplotDecoder
#'
#' A Vplot decoder network.  A DCGAN-like generator was used here.
#'
#' @param vplot_width V-plot width (genomic position-wise)
#' @param vplot_height V-plot height (fragment size wise)
#' @param filters0 The dimensionality of the output space of the first image from the latent space (default: 64L)
#' @param filters The dimensionality of the output space of the deconved images (default: c(32L, 32L, 1L))
#' @param kernel_size  kernel size
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
VplotDecoder <- function(
	vplot_width,	# width of the image
	vplot_height, # height of the image
	filters0 = 64L,
	filters = 32L,
	kernel_size = 3L,
	upsample_layers = 4L,
	strides = c(2L, 2L),
	momentum = 0.8,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vplot_width <- vplot_width
		self$vplot_height <- vplot_height
		self$upsample_layers <- upsample_layers

		stopifnot(vplot_height %% strides[1]^upsample_layers == 0)
		stopifnot(vplot_width %% strides[2]^upsample_layers == 0)

		interval_dim0 <- as.integer(vplot_height / strides[1]^upsample_layers)
		window_dim0 <- as.integer(vplot_width / strides[2]^upsample_layers)
		output_dim0 <- as.integer(window_dim0 * interval_dim0 * filters0)

		self$dense_1 <- tf$keras$layers$Dense(
			units = output_dim0,
			activation = 'relu'
		)

		self$reshape_1 <- tf$keras$layers$Reshape(target_shape = c(interval_dim0, window_dim0, filters0))

		self$deconv <- lapply(1:self$upsample_layers, function(i) tf$keras$Sequential(list(
			tf$keras$layers$Conv2DTranspose(
				filters = filters,
				kernel_size = kernel_size,
				strides = shape(strides[1], strides[2]),
				padding = 'same'
			),
			tf$keras$layers$Activation('relu'),
			tf$keras$layers$BatchNormalization(momentum = momentum)
		)))

		self$conv_final <- tf$keras$layers$Conv2D(
			filters = 1L,
			kernel_size = 1L,
			padding = 'same'
		)

		function(x, ..., training = TRUE){

			x <- x %>%
				self$dense_1() %>%
				self$reshape_1() 

			for (i in 1:self$upsample_layers){
				x <- x %>% 
					self$deconv[[i - 1]]() 
			}

			x <- x %>% 
				self$conv_final()

			x	
		}
	})
}
