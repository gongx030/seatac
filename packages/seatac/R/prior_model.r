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

#' primor_model_gmm
#'
prior_model_gmm <- function(latent_dim, n_components = 10, name = NULL){
	keras_model_custom(name = name, function(self) {

		self$loc <- tf$Variable(matrix(0, n_components, latent_dim), dtype = tf$float32)
		self$raw_scale_diag <- tf$Variable(matrix(0, n_components, latent_dim), dtype = tf$float32)
		self$mixture_logits <- tf$Variable(rep(0, n_components), dtype = tf$float32)

		function (x, mask = NULL) {
			tfd_mixture_same_family(
				components_distribution = tfd_multivariate_normal_diag(loc = self$loc, scale_diag = tf$nn$softplus(self$raw_scale_diag)),
				mixture_distribution = tfd_categorical(logits = self$mixture_logits)
			)
		}
	})
}
