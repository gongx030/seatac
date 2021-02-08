#' VplotDecoder
#'
#' A Vplot decoder network
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
	filters = c(32L, 32L, 1L),
	kernel_size = c(3L, 3L, 3L),
	interval_strides = c(2L, 2L, 1L),
	window_strides = c(2L, 2L, 2L),
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vplot_width <- vplot_width
		self$vplot_height <- vplot_height
		self$filters <- filters
		self$n_layers <- length(filters)

		stopifnot(vplot_width %% prod(window_strides) == 0)
		stopifnot(vplot_height %% prod(interval_strides) == 0)
		stopifnot(self$n_layers == length(kernel_size))
		stopifnot(self$n_layers == length(window_strides))
		stopifnot(self$n_layers == length(interval_strides))

		window_dim0 <- as.integer(vplot_width / prod(window_strides))
		interval_dim0 <- as.integer(vplot_height / prod(interval_strides))
		output_dim0 <- as.integer(window_dim0 * interval_dim0 * filters0)

		self$dense_1 <- tf$keras$layers$Dense(
			units = output_dim0,
			activation = 'relu'
		)

		self$dropout_1 <- tf$keras$layers$Dropout(rate)

		self$reshape_1 <- tf$keras$layers$Reshape(target_shape = c(interval_dim0, window_dim0, filters0))

		self$deconv <- lapply(1:self$n_layers, function(i) 
			if (i == self$n_layers){													
				tf$keras$layers$Conv2DTranspose(
					filters = filters[i],
					kernel_size = kernel_size[i],
					strides = shape(interval_strides[i], window_strides[i]),
					padding = 'same'
				)
			}else{
				tf$keras$layers$Conv2DTranspose(
					filters = filters[i],
					kernel_size = kernel_size[i],
					strides = shape(interval_strides[i], window_strides[i]),
					padding = 'same',
					activation = 'relu'
				)
			}
		)

		self$bn <- lapply(1:self$n_layers, function(i) tf$keras$layers$BatchNormalization())

		function(x, ..., training = TRUE){

			x <- x %>%
				self$dense_1() %>%
				self$dropout_1(training = training) %>%
				self$reshape_1() 

			for (i in 1:self$n_layers){
				x <- x %>% 
					self$deconv[[i - 1]]() %>% # zero-based
					self$bn[[i - 1]]()	# zero-based
			}
			x	
		}
	})
}
