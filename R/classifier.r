#' gmm_cvae_classifier_model
#'
gmm_cvae_classifier_model <- function(
	input_dim,
	name = NULL
){
	keras_model_custom(name = name, function(self){
		self$dense_1 <- layer_dense(units = 1L)
		function(x, mask = NULL){
			y <- x %>% self$dense_1()
			tfd_independent(
				tfd_bernoulli(logits = y)
			)
		}
	})
} # gmm_cvae_classifier_model
