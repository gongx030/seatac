#' MaskLayer
#'
#' A patch extraction layer
#' 
#' @param patch_size Patch size
#' @param num_patches Number of patches
#' @param name model name
#'
MaskLayer <- function(
	ratio = 0.8,
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		function(x, training = FALSE, ...){

			batch_size <- x$shape[[1]]
			n_bins_per_window <- x$shape[[2]]
			h <- x$shape[[3]]

			if (training){
				m <- (tf$random$uniform(shape(batch_size, n_bins_per_window, h)) > ratio) %>%
					tf$cast(tf$float32)
				x <- tf$multiply(x, m)
				x <- x %>% scale_vplot()
			}else{
				x
			}
		}
	})
}


