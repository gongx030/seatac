#' VplotEncoder
#'
VplotEncoder <- function(
	filters = c(32L, 32L, 32L),
	kernel_size = c(3L, 3L, 3L),
	window_strides = c(2L, 2L, 2L),
	interval_strides = c(2L, 2L, 2L),
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$filters <- filters
		self$kernel_size <- kernel_size
		self$window_strides <-  window_strides
		self$interval_strides <-  interval_strides

		self$conv1d_1 <- tf$keras$layers$Conv2D(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(interval_strides[1], window_strides[1]),
			activation = 'relu'
		)

		self$conv1d_2 <- tf$keras$layers$Conv2D(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(interval_strides[2], window_strides[2]),
			activation = 'relu'
		)

		self$conv1d_3 <- tf$keras$layers$Conv2D(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(interval_strides[3], window_strides[3]),
			activation = 'relu'
		)

		self$bn_1 <- tf$keras$layers$BatchNormalization()
		self$bn_2 <- tf$keras$layers$BatchNormalization()
		self$bn_3 <- tf$keras$layers$BatchNormalization()

		self$flatten_1 <- tf$keras$layers$Flatten()

		function(x, training = TRUE, mask = NULL){
			y <- x %>%
			self$conv1d_1() %>% 
			self$bn_1() %>%
			self$conv1d_2() %>% 
			self$bn_2() %>%
			self$conv1d_3() %>% 
			self$bn_3() %>%
			self$flatten_1()
			y
		}
	})
}


#' VplotDecoder
#' 
VplotDecoder <- function(
	vplot_width,	# width of the image
	vplot_height, # height of the image
	filters0 = 64L,
	filters = c(32L, 32L, 1L),
	kernel_size = c(3L, 3L, 3L),
	window_strides = c(2L, 2L, 2L),
	interval_strides = c(2L, 2L, 1L),
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vplot_width <- vplot_width
		self$vplot_height <- vplot_height

		if (vplot_width %% prod(window_strides) != 0)
			stop(sprintf('vplot_width must be a multiple of %d', prod(window_strides)))

		window_dim0 <- as.integer(vplot_width / prod(window_strides))
		interval_dim0 <- as.integer(vplot_height / prod(interval_strides))
		output_dim0 <- as.integer(window_dim0 * interval_dim0 * filters0)

		self$dense_1 <- tf$keras$layers$Dense(
			units = output_dim0,
			activation = 'relu'
		)

		self$dropout_1 <- tf$keras$layers$Dropout(0.8)

		self$deconv_1 <- tf$keras$layers$Conv2DTranspose(
			filters = filters[1],
			kernel_size = kernel_size[1],
			strides = shape(interval_strides[1], window_strides[1]),
			padding = 'same',
			activation = 'relu'
		)

		self$deconv_2 <- tf$keras$layers$Conv2DTranspose(
			filters = filters[2],
			kernel_size = kernel_size[2],
			strides = shape(interval_strides[2], window_strides[2]),
			padding = 'same',
			activation = 'relu'
		)

		self$deconv_3 <- tf$keras$layers$Conv2DTranspose(
			filters = filters[3],
			kernel_size = kernel_size[3],
			strides = shape(interval_strides[3], window_strides[3]),
			padding = 'same'
		)

		self$bn_1 <- tf$keras$layers$BatchNormalization()
		self$bn_2 <- tf$keras$layers$BatchNormalization()

		self$reshape_1 <- tf$keras$layers$Reshape(target_shape = c(interval_dim0, window_dim0, filters0))

		function(x, training = TRUE, mask = NULL){
			y <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3() 
			y
		}
	})
}


#' VaeEncoder
#' 
VaeEncoder <- function(
	latent_dim,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$latent_dim <- latent_dim
		self$vplot_encoder <- VplotEncoder()
		self$dense_1 <- tf$keras$layers$Dense(units = 2 * self$latent_dim)

		function(x, training = TRUE, mask = NULL){

			y <- x %>%
				self$vplot_encoder() %>%
				self$dense_1()

			tfp$distributions$MultivariateNormalDiag(
				loc = y[, 1:self$latent_dim],
				scale_diag = tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3),
				name = 'posterior'
			)
		}
	})
}

#' Decoder
#' 
VaeDecoder <- function(
	vplot_width,	# width of the image
	vplot_height, # height of the image
	filters0 = 64,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vplot_decoder <- VplotDecoder(
			filters0 = filters0,																								 
			vplot_width = vplot_width, 
			vplot_height = vplot_height
		)
	
		function(x, training = TRUE){

			y <- x %>%
				self$vplot_decoder()
	
			x_pred <- y %>% tf$keras$activations$softmax(1L)
			x_pred
		}
	})
}

#' VaeModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
VaeModel <- function(
	 latent_dim = 10L,
	 n_intervals,
	 bin_size = 5L,
	 block_size,
	 filters0 = 128L,
	 fragment_size_range  = c(80L, 320L),
	 fragment_size_interval = 5L,
	 rate = 0.1,

	 name = NULL
){


	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$block_size <- block_size
		self$bin_size <- bin_size
		self$n_bins_per_block <- as.integer(block_size / bin_size)

		self$fragment_size_range <- fragment_size_range
		self$fragment_size_interval <- fragment_size_interval
		br <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		self$breaks <- tf$constant(br)
		self$centers <- tf$constant((br[-1] + br[-length(br)]) / 2)
		
		self$is_nfr <- self$centers <= 100
		self$is_nucleosome <- self$centers >= 180 & self$centers <= 247

		self$encoder <- VaeEncoder(
			latent_dim = latent_dim
		)

		self$decoder <- VaeDecoder(
			filters0 = filters0,
			vplot_width = self$n_bins_per_block,
			vplot_height = n_intervals
		)

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(latent_dim)),
			scale_identity_multiplier = 1
		)

		function(x, training = TRUE){
			posterior <- x %>% self$encoder()
			z <- posterior$sample()
			x_pred <- z %>% self$decoder()

			list(
				posterior = posterior, 
				z = z, 
				vplots = x_pred
			)
		}
	})
}


#' prepare_data
#'
#' Prepare tfdataset for training and testing a VaeModel model
#'
#' @export
#'
setMethod(
	'prepare_data',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		weighted = TRUE
	){

		d <- list()

		y <- assays(x)$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32)

		d$vplots <- scale_vplot(y)

		# add weight for each genomic bin
		w <- tf$reduce_sum(d$vplots, 1L, keepdims = TRUE) > 0
		w <- w %>% tf$cast(tf$float32)
		if (!weighted)
			w <- tf$ones_like(w)

		d$weight <- w

		d <- d %>%
			tensor_slices_dataset()
		d
	}
)

#' fit
#'
#' @export
#'
setMethod(
	'fit',
	signature(
		model = 'VaeModel',
		x = 'tf_dataset'
	),
	function(
		 model,
		 x,
		 batch_size = 32L,
		 epochs = 100L,
		 test_size = 0.15,
		 learning_rate = 1e-4
	 ){

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		train_step <- function(x, w){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x)
				loss_reconstruction <- (tf$squeeze(w, 3L) * reconstrution_loss(x, res$vplots)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction + loss_kl
	 		})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss,
				loss_reconstruction = loss_reconstruction,
				loss_kl = loss_kl
			)
		}

		test_step <- function(x, w){
			res <- model@model(x)
			metric_test <- (tf$squeeze(w, 3L) * reconstrution_loss(x, res$vplots)) %>%
				tf$reduce_sum(shape(1L, 2L)) %>%
				tf$reduce_mean()
			list(
				metric_test = metric_test
			)
		}

		train_step <- tf_function(train_step) # convert to graph mode
		test_step <- tf_function(test_step) # convert to graph mode

		for (epoch in seq_len(epochs)){
			# training 
			loss_train <- NULL 
			loss_train_reconstruction <- NULL
			loss_train_kl <- NULL
			iter <- make_iterator_one_shot(x$train)
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$vplots, batch$weight)
				loss_train <- c(loss_train, as.numeric(res$loss))
				loss_train_reconstruction <- c(loss_train_reconstruction, as.numeric(res$loss_reconstruction))
				loss_train_kl <- c(loss_train_kl, as.numeric(res$loss_kl))
			})

			# testing
			metric_test <- NULL
			iter <- make_iterator_one_shot(x$test)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- test_step(batch$vplots, batch$weight)
				metric_test <- c(metric_test, as.numeric(res$metric_test))
			})

			message(sprintf('epoch=%6.d/%6.d | train_recon_loss=%13.7f | train_kl_loss=%13.7f | train_loss=%13.7f | test_recon_loss=%13.7f', epoch, epochs, mean(loss_train_reconstruction), mean(loss_train_kl), mean(loss_train), mean(metric_test)))

		}
		model
	}
)

#'
setMethod(
	'predict',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 1L, # v-plot per batch
		scale = -10,
		offset = -0.95
	){

		y <- assays(x)$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32)

		res <- model %>% predict(y, batch_size = batch_size, scale = scale, offset = offset)

		x@assays@data$predicted_counts <- res$vplots %>%
			tf$reshape(c(res$vplots$shape[[1]], -1L)) %>%
			as.matrix()
		
		SummarizedExperiment::rowData(x)$latent <- res$z %>% as.matrix()
		SummarizedExperiment::rowData(x)$nucleosome <- res$nucleosome %>% as.matrix()

		x
	}
) # predict
#'



#' 
#' @export
#'
setMethod(
	'predict',
	signature(
		model = 'VaeModel',
		x = 'tensorflow.tensor'
	),
	function(
		model,
		x,
		batch_size = 128L,
		scale = -10,
		offset = -0.95
	){

		z <- list()
		vplots <- list()	
		nucleosome <- list()	

		x <- scale_vplot(x)

		batches <- cut_data(x$shape[[1]], batch_size)

		for (i in 1:length(batches)){
			b <- batches[[i]]
			posterior <- model@model$encoder(x[b, ,  , , drop = FALSE])
			z[[i]] <- posterior$mean()
			x_pred <- model@model$decoder(z[[i]])
			vplots[[i]] <- x_pred
			nucleosome[[i]] <- x_pred %>% vplot2nucleosome(model@model$is_nucleosome, model@model$is_nfr, scale, offset)

		}

		z <- tf$concat(z, axis = 0L)
		vplots <- tf$concat(vplots, axis = 0L)
		nucleosome <- tf$concat(nucleosome, axis = 0L)

		list(vplots = vplots, z = z, nucleosome = nucleosome)
	}
)

#' 
#' @export
#'
setMethod(
	'encode',
	signature(
		model = 'VaeModel',
		x = 'tensorflow.tensor'
	),
	function(
		model,
		x,
		batch_size = 128L
	){

		x <- scale_vplot(x)

		batches <- cut_data(x$shape[[1]], batch_size)
		res <- list()

		for (i in 1:length(batches)){
			b <- batches[[i]]
			posterior <- model@model$encoder(x[b, ,  , , drop = FALSE])
			z <- posterior$mean()
			res[[i]] <- z
		}

		z <- tf$concat(res, axis = 0L)
		z
	}
)

#'
setMethod(
	'encode',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L
	){

		y <- assays(x)$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32)

		z <- model %>% encode(y, batch_size = batch_size)
		SummarizedExperiment::rowData(x)$latent <- z %>% as.matrix()
		x
	}
) # encode

#' 
#' @export
#'
setMethod(
	'encode',
	signature(
		model = 'VaeModel',
		x = 'GRanges'
	),
	function(
		model,
		x,
		filename,
		genome,
		batch_size = 128L,
		batch_size_read_vplot = 1024L
	){

		window_size <- width(x[1])
		stopifnot(all(width(x) == window_size))

		block_size <- model@model$block_size
		bin_size <- model@model$bin_size
		n_bins_per_block <- model@model$n_bins_per_block

		stopifnot(window_size >= block_size)
		stopifnot(window_size %% bin_size == 0)

		x <- x %>%
			resize(width = window_size - block_size, fix = 'center') %>%
			slidingWindows(bin_size, bin_size) %>%
		 	unlist()

		fragment_size_range <- c(model@model$fragment_size_range[[0]], model@model$fragment_size_range[[1]])
		n_reads <- x %>%
		  resize(width = block_size, fix = 'center') %>% 
		  count_reads(filename, genome, fragment_size_range = fragment_size_range)

		mcols(x)$n_reads <- n_reads
		x <- x[x$n_reads > 0]

		x <- x %>%
			resize(width = block_size, fix = 'center')

		batches <- cut_data(length(x), batch_size_read_vplot)

		latent <- matrix(NA, length(x), model@model$encoder$latent_dim)
		for (i in 1:length(batches)){
			message(sprintf('encode | batch=%6.d/%6.d', i, length(batches)))
			b <- batches[[i]]
			xb <- x[b] %>% 
				read_vplot(filename, genome, bin_size, fragment_size_range) 
			xb <- model %>% encode(xb, batch_size)
			latent[b, ] <- rowData(xb)$latent
		}
		mcols(x)$latent <- latent
		x
	}
)


#' 
#' @export
#'
setMethod(
	'decode',
	signature(
		model = 'VaeModel',
		x = 'tensorflow.tensor'
	),
	function(
		model,
		x,
		batch_size = 128L,
		scale = -10,
		offset = -0.95
	){

		batches <- cut_data(x$shape[[1]], batch_size)
		vplots <- list()
		nucleosome <- list()

		for (i in 1:length(batches)){
			b <- batches[[i]]
			x_pred <- model@model$decoder(x[b, , drop = FALSE])
			vplots[[i]] <- x_pred
			nucleosome[[i]] <- x_pred %>% vplot2nucleosome(model@model$is_nucleosome, model@model$is_nfr, scale, offset)
		}

		vplots  <- tf$concat(vplots, axis = 0L)
		nucleosome <- tf$concat(nucleosome, axis = 0L)

		list(vplots = vplots, nucleosome = nucleosome)
	}
)

#'
#' @export 
#'
setMethod(
	'decode',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 1L, # v-plot per batch
		scale = -10,
		offset = -0.95
	){

		if (is.null(rowData(x)$latent))
			stop('latent must be defined')

		z <- rowData(x)$latent %>% 
			tf$cast(tf$float32) 
		res <- decode(model, z, batch_size = batch_size, scale = scale, offset = offset)

		x@assays@data$predicted_counts <- res$vplots %>%
			tf$reshape(c(res$vplots$shape[[1]], -1L)) %>%
			as.matrix()
		
		SummarizedExperiment::rowData(x)$nucleosome <- res$nucleosome %>% as.matrix()
		x
	}
)
