#' regularizer_model_kmeans
#'
regularizer_model_kmeans <- function(
	latent_dim,
	n_components,
	sigma = 1,
	name = NULL
){

  keras_model_custom(name = name, function(self){

		self$centers <- tf$Variable(tf$random$normal(shape(n_components, latent_dim)))

		function(x, mask = NULL){

			u <- tf$expand_dims(x[[1]], 1L)
			c <- tf$expand_dims(self$centers, 0L)
			d <- tf$subtract(u, c) %>% 
				tf$reduce_sum(2L) %>% 
				tf$square()

			tf$reduce_sum(tf$multiply(x[[2]], d))
		}
	})
} # regularizer_model_kmeans

