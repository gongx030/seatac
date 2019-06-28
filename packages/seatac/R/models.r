#' build_model
#'

gmm_vae <- function(input_dim, feature_dim, latent_dim, n_components, num_samples){

	learnable_prior_model <- function(latent_dim, n_components, name = NULL){
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
	latent_prior_model <- learnable_prior_model(latent_dim, n_components)

	encoder_model <- function(latent_dim, name = NULL){

		keras_model_custom(name = name, function(self){
			self$conv_1 <- layer_conv_2d(
				filters = 4L,
				kernel_size = 3L,
				strides = 2L,
				activation = 'relu'
			)
			self$flatten_1 <- layer_flatten()
			self$dense_1 <- layer_dense(units = 2 * latent_dim)

			function(x, mask = NULL){

				y <- x %>% 
					self$conv_1() %>%
					self$flatten_1() %>%
					self$dense_1()

				tfd_multivariate_normal_diag(
					loc = y[, 1:latent_dim],
					scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
				)
			}
		})
  }
  encoder <- encoder_model(latent_dim)

	decoder_model <- function(input_dim, feature_dim, num_samples, name = NULL){

		keras_model_custom(name = name, function(self){

			self$embedding_1 <- layer_embedding(input_dim = num_samples, output_dim = latent_dim, input_length = 1L)
			self$dense_1 <- layer_dense(units = 8 * 8 * 16, activation = "relu")
			self$reshape_1 <- layer_reshape(target_shape = c(8L, 8L, 16L))
			self$deconv_1 <- layer_conv_2d_transpose(
				filters = 16L,
				kernel_size = 3L,
				strides = 4L,
				padding = 'same',
				activation = 'relu'
			)
			self$deconv_3 <- layer_conv_2d_transpose(
				filters = 1L,
				kernel_size = 3L,
				strides = 1L,
				padding = 'same'
			)

			function(x, mask = NULL){

				y <- list(x[[1]], x[[2]] %>% self$embedding_1()) %>%
					layer_concatenate() %>%	
					self$dense_1() %>%
					self$reshape_1() %>%
					self$deconv_1() %>%
					self$deconv_3()
				
				tfd_independent(
					tfd_bernoulli(logits = y), 
					reinterpreted_batch_ndims = 3L
				)
			}
		})
	}
	decoder <- decoder_model(input_dim, feature_dim, num_samples)

	structure(list(
		encoder = encoder, 
		decoder = decoder, 
		latent_prior_model = latent_prior_model, 
		trainable_prior = TRUE,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		n_components = n_components,
		num_samples = num_samples
	), class = c('seatac_model', 'gmm_vae'))

} # gmm_vae


hmm_vae <- function(input_dim, feature_dim, latent_dim, n_components, num_samples){

	learnable_prior_model <- function(latent_dim, n_components, name = NULL){
		keras_model_custom(name = name, function(self) {

#			self$transition_distribution <- tf$get_variable(

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
	latent_prior_model <- learnable_prior_model(latent_dim, n_components)

} # hmm_vae

saveModel <- function(x, dir){

	if (missing(dir))
		stop('dir must be specified')

	if (!file.exists(dir))
		dir.create(dir, recursive = TRUE)

	encoder_file <- sprintf('%s/encoder.h5', dir)
	flog.trace(sprintf('writing %s', encoder_file))
	save_model_weights_hdf5(x$encoder, encoder_file)

	decoder_file <- sprintf('%s/decoder.h5', dir)
	flog.trace(sprintf('writing %s', decoder_file))
	save_model_weights_hdf5(x$decoder, decoder_file)

	latent_prior_model_file <- sprintf('%s/latent_prior_model.h5', dir)
	flog.trace(sprintf('writing %s', latent_prior_model_file))
	save_model_weights_hdf5(x$latent_prior_model, latent_prior_model_file)

	x$encoder_file <- encoder_file
	x$decoder_file <- decoder_file
	x$latent_prior_model_file <- latent_prior_model_file

	x$encoder <- NULL
	x$decoder <- NULL
	x$latent_prior_model <- NULL

	model_file <- sprintf('%s/model.rds', dir)
	flog.trace(sprintf('writing %s', model_file))
	saveRDS(x, model_file)

} # saveModel


loadModel <- function(dir){

	if (missing(dir))
		stop('dir must be specified')

	encoder_file <- sprintf('%s/encoder.h5', dir)
	decoder_file <- sprintf('%s/decoder.h5', dir)
	latent_prior_model_file <- sprintf('%s/latent_prior_model.h5', dir)
	model_file <- sprintf('%s/model.rds', dir)

	if (!file.exists(encoder_file))
		stop(sprintf('%s does not exist', encoder_file))

	if (!file.exists(decoder_file))
		stop(sprintf('%s does not exist', decoder_file))

	if (!file.exists(latent_prior_model_file))
		stop(sprintf('%s does not exist', latent_prior_model_file))

	if (!file.exists(model_file))
		stop(sprintf('%s does not exist', model_file))

	x <- readRDS(model_file)
	model <- gmm_vae(input_dim = x$input_dim, feature_dim = x$feature_dim, latent_dim = x$latent_dim, n_components = x$n_components)

	model$encoder$set_weights(encoder_file)
	model$decoder$set_weights(decoder_file)
	model$latent_prior_model$set_weights(latent_prior_model_file)

	for (y in names(x)[!names(x) %in% c('encoder', 'decoder', 'latent_prior_model')])
		model[[y]] <- x[[y]]

	class(model) <- 'seatac_model'

	model
} # loadModel


