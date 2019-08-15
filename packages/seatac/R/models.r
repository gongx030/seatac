#' vae
#'

vae <- function(input_dim, feature_dim, latent_dim, n_components, num_samples, prior = 'gmm'){

	if (prior == 'gmm'){
		latent_prior_model <- gmm_prior_model(2 * latent_dim + 1, n_components)
	}else
		stop(sprintf('unknown prior model: %s', prior))

	structure(list(
		encoder = list(
			vplot = vplot_encoder_model(latent_dim),
			coverage = coverage_encoder_model(latent_dim)
		),
		classifier = classifier_model(input_dim, feature_dim),
		decoder = list(
			vplot = vplot_decoder_model(input_dim, feature_dim), 
			coverage = coverage_decoder_model(input_dim)
		),
		latent_prior_model = latent_prior_model,
		trainable_prior = TRUE,
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		n_components = n_components,
		num_samples = num_samples,
		prior = prior
	), class = c('seatac_model', 'vae'))

} # vae


saveModel <- function(x, dir){

	if (missing(dir))
		stop('dir must be specified')

	if (!file.exists(dir))
		dir.create(dir, recursive = TRUE)

	encoder_vplot_file <- sprintf('%s/encoder_vplot.h5', dir)
	flog.info(sprintf('writing %s', encoder_vplot_file))
	x$encoder$vplot$save_weights(encoder_vplot_file)
	x$encoder$vplot_weight_file <- encoder_vplot_file

	encoder_coverage_file <- sprintf('%s/encoder_coverage.h5', dir)
	flog.info(sprintf('writing %s', encoder_coverage_file))
	x$encoder$coverage$save_weights(encoder_coverage_file)
	x$encoder$coverage_weight_file <- encoder_coverage_file

	decoder_vplot_file <- sprintf('%s/decoder_vplot.h5', dir)
	flog.info(sprintf('writing %s', decoder_vplot_file))
	x$decoder$vplot$save_weights(decoder_vplot_file)
	x$decoder$vplot_weight_file <- decoder_vplot_file

	decoder_coverage_file <- sprintf('%s/decoder_coverage.h5', dir)
	flog.info(sprintf('writing %s', decoder_coverage_file))
	x$decoder$coverage$save_weights(decoder_coverage_file)
	x$decoder$coverage_weight_file <- decoder_coverage_file

	latent_prior_file <- sprintf('%s/latent_prior_model.h5', dir)
	flog.info(sprintf('writing %s', latent_prior_file))
	x$latent_prior_model$save_weights(latent_prior_file)
	x$latent_prior_file_weight_file <- latent_prior_file

	x$encoder$vplot <- NULL
	x$encoder$coverage <- NULL
	x$decoder$vplot <- NULL
	x$decoder$coverage <- NULL
	x$latent_prior_model <- NULL

	model_file <- sprintf('%s/model.rds', dir)
	flog.info(sprintf('writing %s', model_file))
	saveRDS(x, model_file)

} # saveModel


loadModel <- function(dir){

	if (missing(dir))
		stop('dir must be specified')

	model_file <- sprintf('%s/model.rds', dir)

	if (!file.exists(model_file))
		stop(sprintf('%s does not exist', model_file))

	x <- readRDS(model_file)

	model <- vae(
		input_dim = x$input_dim,
		feature_dim = x$feature_dim,
		latent_dim = x$latent_dim,
		n_components = x$n_components,
		num_samples = x$num_samples,
		prior = x$prior
	)

  # reactivate the layers
  array(0, dim = c(1, model$input_dim, model$feature_dim, 1)) %>%
	  tf$cast(tf$float32) %>%
	  model$encoder$vplot()
	model$encoder$vplot$load_weights(x$encoder$vplot_weight_file)

  array(0, dim = c(1, model$input_dim, 1)) %>%
	  tf$cast(tf$float32) %>%
	  model$encoder$coverage()
	model$encoder$coverage$load_weights(x$encoder$coverage_weight_file)

  array(0, dim = c(1, model$latent_dim * 2)) %>%
	  tf$cast(tf$float32) %>%
	  model$decoder$vplot()
	model$decoder$vplot$load_weights(x$decoder$vplot_weight_file)

  array(0, dim = c(1, model$latent_dim * 2)) %>%
	  tf$cast(tf$float32) %>%
	  model$decoder$coverage()
	model$decoder$coverage$load_weights(x$decoder$coverage_weight_file)

	model$latent_prior_model$load_weights(x$latent_prior_file_weight_file)

	model
} # loadModel


#' encoding model for the vplot
#'
vplot_encoder_model <- function(
	latent_dim, 
	filters = c(32L, 32L, 32L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$conv_1 <- layer_conv_2d(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(input_strides[1], feature_strides[1]),
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$conv_2 <- layer_conv_2d(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			activation = 'relu'
		)

		self$bn_2 <- layer_batch_normalization()

		self$conv_3 <- layer_conv_2d(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
			activation = 'relu'
		)

		self$bn_3 <- layer_batch_normalization()

		self$flatten_1 <- layer_flatten()
		self$dense_1 <- layer_dense(units = 2 * latent_dim)

		function(x, mask = NULL){

			y <- x %>% 
				self$conv_1() %>%
				self$bn_1() %>%
				self$conv_2() %>%
				self$bn_2() %>%
				self$conv_3() %>%
				self$bn_3() %>%
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


#' vplot_decoder_model
#'
#' @param input_dim input dimension of the v-plot (i.e. the width of the genomic region)
#' @param feature_dim feature dimension of the v-plot	(i.e. the fragment size dimension)
#' @param filters0 the beginning filter dimension coming out of the latent layer
#' @param filters the filter sizes of each deconv layer
#' @param kernel_size the kernel size of each deconv layer.  The feature and input spaces shared the same kernel size. 
#'
vplot_decoder_model <- function(
	input_dim, 
	feature_dim, 
	filters0 = 32L, 
	filters = c(32L, 32L, 1L), 
	kernel_size = c(3L, 3L, 3L), 
	input_strides = c(2L, 2L, 2L), 
	feature_strides = c(2L, 2L, 1L), 
	name = NULL
){

	input_dim0 <- input_dim / prod(input_strides)
	feature_dim0 <- feature_dim / prod(feature_strides)
	output_dim0 <- input_dim0 * feature_dim0 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(input_dim0, feature_dim0, filters0))

		self$deconv_1 <- layer_conv_2d_transpose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(input_strides[1], feature_strides[1]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_1 <- layer_batch_normalization()

		self$deconv_2 <- layer_conv_2d_transpose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(input_strides[2], feature_strides[2]),
			padding = 'same',
			activation = 'relu'
		)
		self$bn_2 <- layer_batch_normalization()

		self$deconv_3 <- layer_conv_2d_transpose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(input_strides[3], feature_strides[3]),
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
				self$bn_2() %>%
				self$deconv_3()
			
			tfd_independent(
				tfd_bernoulli(logits = y), 
				reinterpreted_batch_ndims = 3L
			)
		}
	})
}

coverage_decoder_model <- function(input_dim, filters0 = 8L, filters = c(8L, 1L), kernel_size = c(3L, 3L), strides = c(4L, 1L), name = NULL){

	input_dim0 <- input_dim / prod(strides)
	output_dim0 <- input_dim0 * 1 * filters0

	keras_model_custom(name = name, function(self){

		self$dense_1 <- layer_dense(units = output_dim0, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$reshape_1 <- layer_reshape(target_shape = c(input_dim0, 1L, filters0))
		self$reshape_2 <- layer_reshape(target_shape = c(input_dim, 1L))

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

classifier_model <- function(input_dim, feature_dim, name = NULL){

	keras_model_custom(name = name, function(self){

		self$conv_vplot_1 <- layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = shape(2L, 2L)
		)
		self$bn_vplot_1 <- layer_batch_normalization()
		self$activation_vplot_1 <- layer_activation(activation = 'relu')
		self$pooling_vplot_1 <- layer_max_pooling_2d(pool_size = c(2L, 2L))

		self$conv_vplot_2 <- layer_conv_2d(
			filters = 32L,
			kernel_size = 3L,
			strides = shape(2L, 2L)
		)
		self$bn_vplot_2 <- layer_batch_normalization()
		self$activation_vplot_2 <- layer_activation(activation = 'relu')
		self$pooling_vplot_2 <- layer_max_pooling_2d(pool_size = c(2L, 2L))

		self$flatten_vplot_1 <- layer_flatten()

		self$conv_coverage_1 <- layer_conv_1d(
			filters = 8L,
			kernel_size = 3L,
			strides = 2L
		)
		self$bn_coverage_1 <- layer_batch_normalization()
		self$activation_coverage_1 <- layer_activation(activation = 'relu')
		self$pooling_coverage_1 <- layer_max_pooling_1d(pool_size = 4L)

		self$flatten_coverage_1 <- layer_flatten()

		self$dense_1 <- layer_dense(units = 16L, activation = 'relu')
		self$dropout_1 <- layer_dropout(rate = 0.2)
		self$dense_2 <- layer_dense(units = 1L)

		function(x, mask = NULL){

			h_vplot <- x$vplot %>% 

#				self$conv_vplot_1() %>%
#				self$bn_vplot_1() %>%
#				self$activation_vplot_1() %>%
#				self$pooling_vplot_1() %>%

#				self$conv_vplot_2() %>%
#				self$bn_vplot_2() %>%
#				self$activation_vplot_2() %>%
#				self$pooling_vplot_2() %>%

				self$flatten_vplot_1()
			
			h_coverage <- x$coverage %>%

#				self$conv_coverage_1() %>%
#				self$bn_coverage_1() %>%
#				self$activation_coverage_1() %>%
#				self$pooling_coverage_1() %>%

				self$flatten_coverage_1()
			
			y <- layer_concatenate(list(h_vplot, h_coverage)) %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$dense_2()

			tfd_bernoulli(logits = y)
		}
	})
}
