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
	
		self$final_layer <- tf$keras$layers$Dense(
			units = 1L,
			activation = 'sigmoid'
		)
	
		function(x, training = TRUE){

			y <- x %>%
				self$vplot_decoder()
	
			x_pred <- y %>% tf$keras$activations$softmax(1L)
	
			nucleosome <- y %>%
				tf$squeeze(3L) %>%
				tf$transpose(shape(0L, 2L, 1L)) %>%
				self$final_layer() %>%
				tf$squeeze(2L)

			list(x_pred = x_pred, nucleosome = nucleosome)
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
	 rate = 0.1,
	 name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$block_size <- block_size
		self$bin_size <- bin_size
		self$n_bins_per_block <- as.integer(block_size / bin_size)

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
			res <- z %>% self$decoder()

			list(
				posterior = posterior, 
				z = z, 
				x_pred = res$x_pred,
				nucleosome = res$nucleosome
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
		x
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

		w <- y %>% tf$reduce_sum(shape(1L), keepdims = TRUE)  # sum of reads per bin
		y <- y / tf$where(w > 0, w, tf$ones_like(w))  # scale so that the sum of each bar is one (softmax)
		d$vplots <- y

		# add weight for each genomic bin
		w <- tf$reduce_sum(d$vplots, 1L, keepdims = TRUE) > 0
		w <- w %>% tf$cast(tf$float32)
		d$weight <- w

		y <- rowData(x)$nucleoatac %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$expand_dims(2L) %>%
			tf$nn$avg_pool1d(ksize = x@bin_size, strides = x@bin_size, padding = 'VALID') %>%
			tf$squeeze(2L) %>%
			scale01()

		d$nucleoatac  <- y

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
		 test_size = 0.15
	 ){

		optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		train_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')
		bce <- tf$keras$losses$BinaryCrossentropy()

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		train_step <- function(x, w, y){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x)
				loss_reconstruction <- (tf$squeeze(w, 3L) * train_loss(x, res$x_pred)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss_nucleosome <- train_loss(y, res$nucleosome) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction + loss_kl + loss_nucleosome
	 		})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss,
				loss_reconstruction = loss_reconstruction,
				loss_kl = loss_kl,
				loss_nucleosome = loss_nucleosome
			)
		}


		test_step <- function(x, w, y){
			res <- model@model(x)
			metric_nucleoatac <- bce(y, res$nucleosome)
			metric_test <- (tf$squeeze(w, 3L) * train_loss(x, res$x_pred)) %>%
				tf$reduce_sum(shape(1L, 2L)) %>%
				tf$reduce_mean()
			list(
				metric_test = metric_test,
				metric_nucleoatac = metric_nucleoatac
			)
		}

		train_step <- tf_function(train_step) # convert to graph mode
		test_step <- tf_function(test_step) # convert to graph mode

		for (epoch in seq_len(epochs)){
			# training 
			loss_train <- NULL 
			loss_train_reconstruction <- NULL
			loss_train_kl <- NULL
			loss_train_nucleosome <- NULL
			iter <- make_iterator_one_shot(x$train)
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$vplots, batch$weight, batch$nucleoatac)
				loss_train <- c(loss_train, as.numeric(res$loss))
				loss_train_reconstruction <- c(loss_train_reconstruction, as.numeric(res$loss_reconstruction))
				loss_train_kl <- c(loss_train_kl, as.numeric(res$loss_kl))
				loss_train_nucleosome <- c(loss_train_nucleosome, as.numeric(res$loss_nucleosome))
			})

			# testing
			metric_test <- NULL
			metric_nucleoatac <- NULL
			iter <- make_iterator_one_shot(x$test)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- test_step(batch$vplots, batch$weight, batch$nucleoatac)
				metric_test <- c(metric_test, as.numeric(res$metric_test))
				metric_nucleoatac <- c(metric_nucleoatac, as.numeric(res$metric_nucleoatac))
			})

			message(sprintf('epoch=%6.d/%6.d | train_recon_loss=%13.7f | train_kl_loss=%13.7f | train_nucleosome_loss=%13.7f | train_loss=%13.7f | test_recon_loss=%13.7f | test_nucleosome=%13.7f', epoch, epochs, mean(loss_train_reconstruction), mean(loss_train_kl), mean(loss_train_nucleosome), mean(loss_train), mean(metric_test), mean(metric_nucleoatac)))

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
		batch_size = 1L # v-plot per batch
	){

		batches <- cut_data(length(x), batch_size)

		n_blocks_per_window <- as.integer(x@n_bins_per_window - model@model$n_bins_per_block + 1)
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)
		predicted_nucleosome <- matrix(0, length(x), x@n_bins_per_window)

		predict_step <- function(x){
			posterior <- x %>% model@model$encoder()
			res <- posterior$mean() %>% model@model$decoder()
			res
		}

		for (i in 1:length(batches)){
			if (i == 1 || i %% 100 == 0){
				message(sprintf('predicting | batch=%4.d/%4.d', i, length(batches)))
			}

			b <- batches[[i]]
			d <- select_blocks(
				x[b],
				block_size = model@model$block_size,
				min_reads = 0,
				with_kmers = FALSE
			)

			res <- predict_step(d$vplots)

			# whether the block has read
			# the blocks without any reads will be ignored when aggregating the blocks for a window
			# note that VAE will give a V-plot prediction even the block has no reads at all
			include <- tf$cast(d$n > 0, tf$float32)

			w <- tf$reshape(include, shape(d$n$shape[[1]], 1L, 1L, 1L))
			x_pred <- tf$multiply(d$vplots, 1 - w) + tf$multiply(res$x_pred, w)

			x_pred <- x_pred %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, model@model$n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			x_pred <- aperm(x_pred, c(1, 3, 2))
			dim(x_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)
			predicted_counts[b, ] <- x_pred

			nucleosome0 <- tf$zeros_like(res$nucleosome)
			w <- tf$reshape(include, shape(d$n$shape[[1]], 1L))

			nucleosome_pred <- tf$multiply(nucleosome0, 1 - w) + tf$multiply(res$nucleosome, w)

			predicted_nucleosome[b, ] <- nucleosome_pred %>%
				tf$expand_dims(2L) %>%
				tf$expand_dims(3L) %>%
				tf$transpose(shape(0L, 2L, 1L, 3L)) %>% 
				tf$reshape(c(length(b), n_blocks_per_window, 1L, model@model$n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks()  %>%
				tf$squeeze(shape(1L, 3L)) %>%
				as.matrix()
		}

		SummarizedExperiment::assays(x)$predicted_counts <- predicted_counts
		SummarizedExperiment::rowData(x)$predicted_nucleosome <- predicted_nucleosome
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
		batch_size = 128L
	){

		res <- list()
		batches <- cut_data(length(x), batch_size)

		for (i in 1:length(batches)){
			b <- batches[[i]]
			res[[i]] <- model@model(x[b, , , ])
		}

		x_pred <- tf$concat(lapply(res, function(r) r$x_pred), axis = 0L)
		nucleosome <- tf$concat(lapply(res, function(r) r$nucleosome), axis = 0L)
		z <- tf$concat(lapply(res, function(r) r$z), axis = 0L)

		list(vplots = vplots, nucleosome = nucleosome, z = z)
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

		batches <- cut_data(x$shape[[1]], batch_size)
		res <- list()

		for (i in 1:length(batches)){
			b <- batches[[i]]
			posterior <- model@model$encoder(x[b, ,  , ])
			z <- posterior$mean()
			res[[i]] <- z
		}

		z <- tf$concat(res, axis = 0L)
		z
	}
)

