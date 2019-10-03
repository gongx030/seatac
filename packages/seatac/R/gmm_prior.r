#' gmm_prior_model
#'
gmm_prior_model <- function(latent_dim, n_components, name = NULL){
	keras_model_custom(name = name, function(self) {
		self$loc <- tf$get_variable(
			name = 'loc',
			shape = list(n_components, latent_dim),
			dtype = tf$float32
		)
		self$raw_scale_diag <- tf$get_variable(
			name = 'raw_scale_diag',
			shape = list(n_components, latent_dim),
			dtype = tf$float32
		)
		self$mixture_logits <- tf$get_variable(
			name = 'mixture_logits',
			shape = c(n_components),
			dtype = tf$float32
		)
		function (x, mask = NULL) {
			tfd_mixture_same_family(
				components_distribution = tfd_multivariate_normal_diag(loc = self$loc, scale_diag = tf$nn$softplus(self$raw_scale_diag)),
				mixture_distribution = tfd_categorical(logits = self$mixture_logits)
			)
		}
	})
}

