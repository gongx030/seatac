#' prior_model
#'
prior_model <- function(latent_dim, name = NULL){
	keras_model_custom(name = name, function(self) {
		function (x, mask = NULL) {
			tfd_multivariate_normal_diag(
				loc  = tf$zeros(list(latent_dim)),
				scale_identity_multiplier = 1
			)
		}
	})
} # prior model
