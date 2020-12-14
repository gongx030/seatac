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

		function(x, fragment_size, training = TRUE, mask = NULL){

			y <- x %>%
				self$vplot_encoder()

			y <- tf$concat(list(y, fragment_size), axis = 1L)
			y <- y %>%	self$dense_1()

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

		function(x, fragment_size, training = TRUE){

			posterior <- self$encoder(x, fragment_size)
			z <- posterior$sample()

			z_cond <- tf$concat(list(z, fragment_size), axis = 1L)
			x_pred <- z_cond %>% self$decoder()

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
		x
	){

		d <- list()

		d$vplots <- assays(x)$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32) 

		d$vplots <- d$vplots %>%
			scale_vplot()

		# add weight for each genomic bin
		w <- tf$reduce_sum(d$vplots, 1L, keepdims = TRUE) > 0
		w <- w %>% tf$cast(tf$float32)
		d$weight <- w

		d
	}
)


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
		x = 'VplotsList'
	),
	function(
		model,
		x
	){

		d <- lapply(1:length(x), function(i){
			model %>% prepare_data(x[[i]])
		})

		for (j in names(d[[1]])){
			d[[1]][[j]] <- tf$concat(lapply(1:length(x), function(i) d[[i]][[j]]), axis = 0L)
		}
		d[[1]]
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
		 learning_rate = 1e-4,
		 num_reads = NA
	 ){

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		train_step <- function(x, w){

			fragment_size <- get_fragment_size(x)

			if (is.na(num_reads)){
				x_input <- x
			}else{
				x_input <- downsample_vplot(x, num_reads)
			}

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x_input, fragment_size)
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
			fragment_size <- get_fragment_size(x)
			res <- model@model(x, fragment_size)
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
			iter <- x$train %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
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
#' @export 
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
		step = 20L,
		scale = -10,
		offset = -0.95,
		min_reads = 1
	){

		stopifnot(step %% x@bin_size == 0)

		# select the windows that have at least one read
		n <- rowSums(assays(x)$counts)
		x <- x[n > 0]

		block_size <- model@model$block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_bins_per_step <- as.integer(step / x@bin_size)

		block_starts <- seq(1, x@n_bins_per_window - n_bins_per_block + 1, by = n_bins_per_step)
		block_ends <- block_starts + n_bins_per_block - 1
		n_blocks_per_window <- length(block_starts)

		batch_size_window <- floor(batch_size / n_blocks_per_window)
		batches <- cut_data(length(x), batch_size_window)

		nucleosome <- list()
		blocks <- list()
		xb_queue <- tf$zeros(shape(0L,  x@n_intervals, n_bins_per_block, 1L))

		fragment_size <- assays(x)$counts %>% colSums() %>%
			matrix(x@n_bins_per_window, x@n_intervals) %>%
			colSums()
		fragment_size <- fragment_size / sum(fragment_size) 
		fragment_size <- tf$cast(fragment_size, tf$float32)

		pred_step <- function(x){
			fs <- fragment_size %>%
				tf$expand_dims(0L) %>%
				tf$`repeat`(repeats = x$shape[[1]], axis = 0L)
			x <- x %>% scale_vplot() 
			posterior <- model@model$encoder(x, fs)
			z <- posterior$mean()
			z_cond <- tf$concat(list(z, fs), axis = 1L)
			xb_pred <- model@model$decoder(z_cond)
			xb_pred
		}
		pred_step <- tf_function(pred_step) # convert to graph mode

		for (i in 1:length(batches)){

			message(sprintf('predict | batch=%6.d/%6.d', i, length(batches)))

			b <- batches[[i]]

			xb <- assays(x[b])$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					x@n_intervals,
					x@n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
		  	tf$image$extract_patches(
					sizes = c(1L, x@n_intervals, n_bins_per_block, 1L),
					strides = c(1L, 1L, n_bins_per_step, 1L),
					rates = c(1L, 1L, 1L, 1L),
					padding = 'VALID'
				) %>%
				tf$squeeze(axis = 1L) %>%
			 	tf$reshape(c(length(b) * n_blocks_per_window, x@n_intervals, n_bins_per_block)) %>%
				tf$expand_dims(-1L) 

			n_reads <- xb %>% tf$reduce_sum(shape(1L, 2L, 3L))
			valid <- n_reads >= min_reads

			xb <- xb %>% tf$boolean_mask(valid)

			blocks[[i]] <- x[b] %>%
				granges() %>%
				slidingWindows(width = block_size, step = step) %>% 
				unlist()

			blocks[[i]] <- blocks[[i]][as.logical(valid)]

			if (xb$shape[[1]] > 0){
				xb_queue <- tf$concat(list(xb_queue, xb), axis = 0L)
			}

			if (xb_queue$shape[[1]] >= batch_size){

				update <- rep(FALSE, xb_queue$shape[[1]])
				update[1:batch_size] <- TRUE

				xb_pred <- xb_queue %>% tf$boolean_mask(update) %>% pred_step()
				nucleosome <- c(nucleosome, xb_pred %>% vplot2nucleosome(model@model$is_nucleosome, model@model$is_nfr, scale, offset))
				xb_queue <- xb_queue %>% tf$boolean_mask(!update)
			}

		}

		if (xb_queue$shape[[1]] > 0){	# predict the remaining vplots in the queue
			xb_pred <- xb_queue %>% pred_step()
			nucleosome <- c(nucleosome, xb_pred %>% vplot2nucleosome(model@model$is_nucleosome, model@model$is_nfr, scale, offset))
		}

		if (length(nucleosome) > 0){

			nucleosome <- tf$concat(nucleosome, axis = 0L)

			blocks <- blocks %>% 
				GRangesList() %>%
				unlist()

			nucleosome <- nucleosome %>%
				tf$reshape(shape(-1L)) %>%
				as.numeric() 
		
			nucleosome <- round(nucleosome * 100)

			bins <- blocks %>%
				slidingWindows(x@bin_size, x@bin_size) %>%
				unlist()

			mcols(bins)$score <- nucleosome
			
			cvg <- bins %>% coverage(weight = 'score')
			n <- bins %>% coverage()
			y <- as(cvg / n, 'GRanges')
			y$score[is.na(y$score)] <- 0
			seqinfo(y)@genome <- seqinfo(bins)@genome

			y
		}

	}
) # predict
#'
