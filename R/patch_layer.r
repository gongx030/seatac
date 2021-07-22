#' PatchLayer
#'
#' A patch extraction layer
#' 
#' @param patch_size Patch size
#' @param num_patches Number of patches
#' @param name model name
#'
PatchLayer <- function(
	patch_size = NULL,
	step_size = NULL,
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$patch_size <- patch_size
		self$step_size <- step_size

		function(x, ...){

			batch_size <- x$shape[[1]]
			n_intervals <- x$shape[[2]]

			patches <- tf$image$extract_patches(
				images = x,
				sizes = shape(1L, n_intervals, self$patch_size, 1L),
				strides = shape(1L, 1L, self$step_size, 1L),
				rates = shape(1L, 1L, 1L, 1L),
				padding = "VALID"
			) %>%
				tf$squeeze(axis = 1L)

			patches

		}
	})
}


