#' build_model
#'

vae <- function(input_dim, feature_dim, latent_dim, n_components, num_samples, prior = 'gmm'){

	if (prior == 'gmm'){
		latent_prior_model <- gmm_prior_model(latent_dim, n_components)
	}else if (prior == 'hmm'){
		latent_prior_model <- hmm_prior_model(latent_dim, n_components, num_samples)
	}else
		stop(sprintf('unknown prior model: %s', prior))

	structure(list(
		encoder = list(
			vplot = vplot_encoder_model(latent_dim),
			coverage = coverage_encoder_model(latent_dim)
		),
		decoder = list(
			vplot = vplot_decoder_model(input_dim, feature_dim), 
			coverage = coverage_decoder_model(feature_dim)
		),
		latent_prior_model = list(
			vplot = latent_prior_model, 
			coverage = vanilla_prior_model(latent_dim)
		),
		trainable_prior = TRUE,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		n_components = n_components,
		num_samples = num_samples,
		prior = prior
	), class = c('seatac_model', 'vae'))

} # vae


hmm_prior_model <- function(latent_dim, n_components, num_steps){

	initial_state_logits <- rep(0, n_components)

	transition_prob <- matrix(0.1, n_components, n_components)
	diag(transition_prob) <- 0
	diag(transition_prob) <- 1 - rowSums(transition_prob)
	transition_distribution <- tfd_categorical(
		probs = transition_prob %>%
		tf$cast(tf$float32)
	)

	learnable_prior_model <- function(latent_dim, n_components, num_steps, name = NULL){

		keras_model_custom(name = name, function(self) {
											 
			self$logits <- tf$get_variable(
				name = 'logits', 
				shape = list(n_components, n_components),
				dtype = tf$float32
			)

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

			function (x, mask = NULL) {
				tfd_hidden_markov_model(
					initial_distribution = tfd_categorical(logits = initial_state_logits),
					transition_distribution = transition_distribution,
#					transition_distribution = tfd_categorical(logits = self$logits),
					observation_distribution =  tfd_multivariate_normal_diag(loc = self$loc, scale_diag = tf$nn$softplus(self$raw_scale_diag)),
					num_steps = num_steps
				)
			}
		})
	}
	latent_prior_model <- learnable_prior_model(latent_dim, n_components, num_steps)

} # hmm_prior_model


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


vplot_encoder_model <- function(latent_dim, filters = c(4L), kernel_size = c(3L), strides = c(2L), name = NULL){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = strides[1],
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$flatten_1 <- layer_flatten()
		self$dense_1 <- layer_dense(units = 2 * latent_dim)
		self$dropout_1 <- layer_dropout(rate = 0.2)

		function(x, mask = NULL){

			y <- x %>% 
				self$conv_1() %>%
				self$bn_1() %>%
				self$flatten_1() %>%
				self$dense_1()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})
}

#' encoding model for the coverage data
#
coverage_encoder_model <- function(latent_dim, filters = c(8L), kernel_size = c(3L), strides = c(2L), name = NULL){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_1d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = strides[1],
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$flatten_1 <- layer_flatten()
		self$dense_1 <- layer_dense(units = 2 * latent_dim)
		self$dropout_1 <- layer_dropout(rate = 0.2)

		function(x, mask = NULL){

			y <- x %>% 
				self$conv_1() %>%
				self$bn_1() %>%
				self$flatten_1() %>%
				self$dense_1()

			tfd_multivariate_normal_diag(
				loc = y[, 1:latent_dim],
				scale_diag = tf$nn$softplus(y[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
			)
		}
	})
}


vplot_decoder_model <- function(input_dim, feature_dim, filters0 = 8L, filters = c(8L, 1L), kernel_size = c(3L, 3L), strides = c(4L, 1L), name = NULL){

	input_dim0 <- input_dim / prod(strides)
	feature_dim0 <- feature_dim / prod(strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(input_dim0, feature_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = strides[1],
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = strides[2],
			padding = 'same'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2()
			
			tfd_independent(
				tfd_bernoulli(logits = y), 
				reinterpreted_batch_ndims = 3L
			)
		}
	})
}

coverage_decoder_model <- function(feature_dim, filters0 = 8L, filters = c(8L, 1L), kernel_size = c(3L, 3L), strides = c(4L, 1L), name = NULL){

	input_dim0 <- feature_dim/ prod(strides)
	output_dim0 <- input_dim0 * 1 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(input_dim0, 1L, filters0))
		self$reshape_2 <- layer_reshape(target_shape = c(feature_dim, 1L))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = shape(kernel_size[1], 1L),
			strides = shape(strides[1], 1L),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = shape(kernel_size[2], 1L),
			strides = shape(strides[2], 1L),
			padding = 'same'
		)

		function(x, mask = NULL){

			y <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$reshape_2() 
			
			tfd_independent(
				tfd_bernoulli(logits = y),
				reinterpreted_batch_ndims = 2L
			)
		}
	})
}

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

vanilla_prior_model <- function(latent_dim, name = NULL){
	keras_model_custom(name = name, function(self){
		function(x, mask = NULL){
			tfd_multivariate_normal_diag(
				loc = tf$zeros(shape(latent_dim)),
				scale_identity_multiplier = 1
			)
		}
	})	
}


