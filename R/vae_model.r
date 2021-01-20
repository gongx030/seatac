#' VplotEncoder
#'
#' A Vplot encoder network
#' 
#' @param filters filters
#' @param kernel_size kernel size
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param name model name
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
		self$n_layers <- length(filters)

		self$conv1d <- lapply(1:self$n_layers, function(i) 
			tf$keras$layers$Conv2D(
				filters = filters[i],
				kernel_size = kernel_size[i],
				strides = shape(interval_strides[i], window_strides[i]),
				activation = 'relu'
			)
		)

		self$bn <- lapply(1:self$n_layers, function(i) tf$keras$layers$BatchNormalization())

		function(x, training = TRUE, mask = NULL){

			for (i in 1:self$n_layers){
				x <- x %>% 
					self$conv1d[[i - 1]]() %>% # zero-based
					self$bn[[i - 1]]()	# zero-based
			}
			x
		}
	})
}


#' VplotDecoder
#'
#' A Vplot decoder network
#'
#' @param vplot_width V-plot width (genomic position-wise)
#' @param vplot_height V-plot height (fragment size wise)
#' @param filters0 The dimensionality of the output space of the first image from the latent space (default: 64L)
#' @param filters The dimensionality of the output space of the deconved images (default: c(32L, 32L, 1L))
#' @param kernel_size  kernel size
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
VplotDecoder <- function(
	vplot_width,	# width of the image
	vplot_height, # height of the image
	filters0 = 64L,
	filters = c(32L, 32L, 1L),
	kernel_size = c(3L, 3L, 3L),
	interval_strides = c(2L, 2L, 1L),
	window_strides = c(2L, 2L, 2L),
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vplot_width <- vplot_width
		self$vplot_height <- vplot_height
		self$filters <- filters
		self$n_layers <- length(filters)

		stopifnot(vplot_width %% prod(window_strides) == 0)
		stopifnot(vplot_height %% prod(interval_strides) == 0)
		stopifnot(self$n_layers == length(kernel_size))
		stopifnot(self$n_layers == length(window_strides))
		stopifnot(self$n_layers == length(interval_strides))

		window_dim0 <- as.integer(vplot_width / prod(window_strides))
		interval_dim0 <- as.integer(vplot_height / prod(interval_strides))
		output_dim0 <- as.integer(window_dim0 * interval_dim0 * filters0)

		self$dense_1 <- tf$keras$layers$Dense(
			units = output_dim0,
			activation = 'relu'
		)

		self$dropout_1 <- tf$keras$layers$Dropout(rate)

		self$reshape_1 <- tf$keras$layers$Reshape(target_shape = c(interval_dim0, window_dim0, filters0))

		self$deconv <- lapply(1:self$n_layers, function(i) 
			if (i == self$n_layers){													
				tf$keras$layers$Conv2DTranspose(
					filters = filters[i],
					kernel_size = kernel_size[i],
					strides = shape(interval_strides[i], window_strides[i]),
					padding = 'same'
				)
			}else{
				tf$keras$layers$Conv2DTranspose(
					filters = filters[i],
					kernel_size = kernel_size[i],
					strides = shape(interval_strides[i], window_strides[i]),
					padding = 'same',
					activation = 'relu'
				)
			}
		)

		self$bn <- lapply(1:self$n_layers, function(i) tf$keras$layers$BatchNormalization())

		function(x, ...){

			x <- x %>%
				self$dense_1() %>%
				self$dropout_1() %>%
				self$reshape_1() 

			for (i in 1:self$n_layers){
				x <- x %>% 
					self$deconv[[i - 1]]() %>% # zero-based
					self$bn[[i - 1]]()	# zero-based
			}
			x	
		}
	})
}


#' VaeEncoder
#'
#' A VAE encoder network
#' 
#' @param latent_dim Latent dimension (default: 10L)
#' @param filters Filters
#' @param kernel_size Kernel sizes
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param distribution Output distributions (MultivariateNormalDiag, LogNormal, or None)
#' @param rate Dropout rate (default: 0.1)
#' @param name model name
#' 
VaeEncoder <- function(
	latent_dim = 10L,
	filters = c(32L, 32L, 32L),
	kernel_size = c(3L, 3L, 3L),
	window_strides = c(2L, 2L, 2L),
	interval_strides = c(2L, 2L, 2L),
	distribution = 'MultivariateNormalDiag',
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$latent_dim <- latent_dim

		self$vplot_encoder <- VplotEncoder(
			filters = filters,
			kernel_size = kernel_size,
			window_strides = window_strides,
			interval_strides = interval_strides
		)

		self$flatten_1 <- tf$keras$layers$Flatten()

		if (distribution == 'MultivariateNormalDiag')
			self$dense_1 <- tf$keras$layers$Dense(units = 2 * self$latent_dim)
		else if (distribution == 'LogNormal')
			self$dense_1 <- tf$keras$layers$Dense(units = 2 * self$latent_dim)
		else if (distribution == 'None')
			self$dense_1 <- tf$keras$layers$Dense(units = self$latent_dim)
		else
			stop(sprintf('unknown distribution: %s', distribution))

		function(x, training = TRUE, mask = NULL){

			y <- x %>%
				self$vplot_encoder() %>%
				self$flatten_1()

			y <- y %>%self$dense_1()

			if (distribution == 'MultivariateNormalDiag'){
 	     	q_m <- y[, 1:self$latent_dim]
				q_v <- tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3)
				tfp$distributions$MultivariateNormalDiag(
					loc = q_m,
					scale_diag = q_v
				)
			}else if (distribution == 'LogNormal'){
 	     	q_m <- y[, 1:self$latent_dim]
				q_v <- tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3)
				xx <- tfp$distributions$LogNormal(
					loc = q_m,
					scale = sqrt(q_v)
				)
			}else if (distribution == 'None')
				y
				
		}
	})
}

#' VaeDecoder
#' 
#' A VAE decoder network
#' 
#' @param vplot_width Vplot width
#' @param vplot_height Vplot height
#' @param filters0 The dimensionality of the output space of the first image from the latent space (default: 64L)
#' @param filters The dimensionality of the output space of the deconved images (default: c(32L, 32L, 1L))
#' @param kernel_size  kernel size
#' @param window_strides Position level strides
#' @param interval_strides Fragment size level strides
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#' 
VaeDecoder <- function(
	vplot_width,	# width of the image
	vplot_height, # height of the image
	filters0 = 64,
	filters = c(32L, 32L, 1L),
	kernel_size = c(3L, 3L, 3L),
	window_strides = c(2L, 2L, 2L),
	interval_strides = c(2L, 2L, 1L),
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vplot_decoder <- VplotDecoder(
			vplot_width = vplot_width, 
			vplot_height = vplot_height,
			filters0 = filters0,																								 
			filters = filters,
			kernel_size = kernel_size,
			window_strides = window_strides,
			interval_strides = interval_strides
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
#' A VAE model for V-plot
#'
#' @param latent_dim Latent dimension (default: 10L)
#' @param block_size Block size in base pairs (default: 640L)
#' @param bin_size Bin size in base pairs(default: 5L) 
#' @param filters0 Filter size after the latent layer (default: 128L)
#' @param fragment_size_range  Fragment size ranges (default: c(80L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 5L)
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
VaeModel <- function(
	 latent_dim = 10L,
	 block_size = 640L,
	 bin_size = 5L,
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
		self$n_intervals <- length(br) - 1L
		self$breaks <- tf$constant(br)
		self$centers <- tf$constant((br[-1] + br[-length(br)]) / 2)
		self$positions <- tf$cast(seq(0 + bin_size / 2, block_size - bin_size / 2, by = bin_size) - (block_size / 2), tf$float32)
		
		self$encoder <- VaeEncoder(
			latent_dim = latent_dim,
			filters = c(32L, 32L, 32L),
			kernel_size = c(3L, 3L, 3L),
			window_strides = c(2L, 2L, 2L),
			interval_strides = c(2L, 2L, 2L),
			distribution = 'MultivariateNormalDiag',
			rate = rate
		)

		self$decoder <- VaeDecoder(
			vplot_width = self$n_bins_per_block,
			vplot_height = self$n_intervals,
			filters0 = filters0,
			filters = c(32L, 32L, 1L),
			kernel_size = c(3L, 3L, 3L),
			window_strides = c(2L, 2L, 2L),
			interval_strides = c(2L, 2L, 1L),
			rate = rate
		)

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(latent_dim)),
			scale_identity_multiplier = 1
		)

		function(x, b, training = TRUE){

			posterior <- x %>% self$encoder()
			z <- posterior$sample()
			c <- list(z, b) %>% tf$concat(axis = 1L)
			x_pred <- c %>% self$decoder()

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
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a Vplots object
#'
#' @return a list that include `vplots`, `weight` and `batch`
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
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

		batch <- rowData(x)$batch %>%
			factor(x@samples) %>%
			as.numeric() %>% 
			tf$cast(tf$int32) %>%
			tf$one_hot(x@n_samples)

		x <- assays(x)$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
			scale_vplot()

		w <- tf$reduce_sum(x, 1L, keepdims = TRUE) > 0
		w <- w %>% tf$squeeze(3L)
		w <- w %>% tf$cast(tf$float32)

		list(
			vplots = x,
			weight = w,
			batch = batch
		)
	}
)


#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param min_reads Mininum number of reads of the V-plots that are used for training (default: 25L)
#' @param max_training_samples Maximum training samples (default: 10000L)
#' @param step_size Step size for sub-sampling Vplots
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @return a VaeModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'fit',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		 model,
		 x,
		 batch_size = 256L,
		 epochs = 500L,
		 learning_rate = 1e-4,
		 warmup = 50L,
		 min_reads = 25L,	# min reads for the training blocks
		 max_training_samples = 10000L,
		 step_size = 20L,
		 compile = TRUE
	){
		
		stopifnot(x@window_size >= model@model$block_size)

		x <- x %>% 
			slidingWindows(width = model@model$block_size, step = step_size, batch_size = 4096L)

		counts <- rowSums(assays(x)$counts)
		sprintf('peaks >= %d reads: %d/%d', min_reads, sum(counts >= min_reads), length(x)) %>% 
			message()
		x <- x[counts >= min_reads]

		if (length(x) > max_training_samples){
			sprintf('downsampling %d peaks for training', max_training_samples) %>% 
				message()
			x <- sample(x, max_training_samples)
		}

		x <- model %>% prepare_data(x)
		x <- x %>% tensor_slices_dataset()
		model <- model %>% fit(x, batch_size = batch_size, epochs = epochs, warmup = warmup, compile = compile, learning_rate = learning_rate)
		model
	}
)


#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @return a VaeModel
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
		 learning_rate = 1e-4,
		 warmup = 50L,
		 compile = TRUE
	 ){

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		train_step <- function(x, w, b, beta){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, b)
				loss_reconstruction <- (w * reconstrution_loss(x, res$vplots)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction + beta * loss_kl
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

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		beta <- c(seq(0, 1, length.out = warmup), rep(1, epochs - warmup))

		for (epoch in seq_len(epochs)){
			loss <- NULL 
			loss_reconstruction <- NULL
			loss_kl <- NULL
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$vplots, batch$weight, batch$batch, beta[epoch])
				loss <- c(loss, as.numeric(res$loss))
				loss_reconstruction <- c(loss_reconstruction, as.numeric(res$loss_reconstruction))
				loss_kl <- c(loss_kl, as.numeric(res$loss_kl))
			})

			sprintf('epoch=%6.d/%6.d | beta=%9.3f | recon_loss=%13.7f | kl_loss=%13.7f | loss=%13.7f', epoch, epochs, beta[epoch], mean(loss_reconstruction), mean(loss_kl), mean(loss)) %>%
				message()

		}
		model
	}
)

#' predict
#' 
#' Predict nucleosome scoreo
#' 
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#' @param step Step size for sub-sampling Vplots
#' @param scale Scale factor for calculating the nucleosome score (default: -10)
#' @param offset Offset factor for calculating the nucleosome score (default: -0.95)
#' @param min_reads Mininum number of reads of the V-plot that are used for predicting (default: 1L)
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


#' test_accessibility
#'
#' Testing the accessibility changes between two samples
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param contrasts This argument specifies what comparison to extract from the object to build a results table.
#' 				This argument must be a character vector with exactly three elements: the name of a factor in the design formula, 
#'				the name of the numerator level for the fold change, and the name of the denominator level for the fold change 
#' @param min_reads The mininum number of reads of the Vplots that are used for comparison (default 5L).  
#'				Both Vplots must have at least this amount of reads
#' @param center_size The center region of the Vplots considered for testing (default: 100L)
#' @param sampling Number of sampling used to computing Bayes factor (default: 200L)
#' @param batch_size Batch size (default: 4L)
#'
#' @return a Vplots object that includes the testing results in rowData fields
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
#' 
setMethod(
	'test_accessibility',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		contrasts = NULL,
		min_reads = 5L,
		center_size = 100L,
		sampling = 100L,
		batch_size = 4L
	){

		stopifnot(is.character(contrasts) && length(contrasts) == 3)

		field <- contrasts[1]
		stopifnot(!is.null(mcols(x)[[field]]) && all(mcols(x)[[field]] %in% x@samples))

		stopifnot(all(contrasts[2:3] %in% x@samples && contrasts[2] != contrasts[3]))

		group1 <- mcols(x)[[field]] %in% contrasts[2]
		group2 <- mcols(x)[[field]] %in% contrasts[3]

		stopifnot(all(granges(x[group1]) == granges(x[group2])))

		is_center <- model@model$positions >= -center_size / 2 & model@model$positions <= center_size / 2
    is_nfr <- model@model$centers <= 100
    is_nucleosome <- model@model$centers >= 180 & model@model$centers <= 247

		counts <- matrix(rowSums(assays(x)$counts), nrow = length(x) / x@n_samples, ncol = x@n_samples, dimnames = list(NULL, x@samples))
		include <- rowSums(counts[, contrasts[2:3]] >= min_reads) == 2L
		n <- sum(include)
		include <- rep(include, x@n_samples)

		y <- assays(x[include])$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window)) %>%
			tf$expand_dims(3L) %>%
			scale_vplot() %>%
			tf$reshape(shape(2L, n, x@n_intervals, x@n_bins_per_window, 1L))	

		batch <- rowData(x[include])$batch %>%
			factor(x@samples) %>%
			as.numeric() %>% 
			tf$cast(tf$int32) %>%
			tf$one_hot(x@n_samples) %>%
			tf$reshape(shape(2L, n, -1L))

		bs <- cut_data(n, batch_size)
		S <- list()
		for (j in 1:length(bs)){

			b <- bs[[j]]

			posterior <- y[, b, , , , drop = FALSE] %>%
				tf$reshape(shape(2L * length(b), x@n_intervals, x@n_bins_per_window, 1L)) %>%
				model@model$encoder()

			z <- posterior$sample(sampling)

			bb <- batch[, b, , drop = FALSE] %>% 
				tf$reshape(shape(2L * length(b), -1L)) %>% 
				tf$expand_dims(0L) %>%
				tf$tile(shape(sampling, 1L, 1L))

			si <- list(z, bb) %>% 
				tf$concat(axis = 2L) %>% 
				tf$reshape(shape(sampling* 2L * length(b), -1L)) %>%
				model@model$decoder() %>% 
				vplot2nucleosome(is_nucleosome, is_nfr) %>% 
				tf$boolean_mask(is_center, axis = 1L) %>% 
				tf$reduce_sum(1L, keepdims = TRUE) %>%
				tf$reshape(shape(sampling, 2L, length(b)))

			S[[j]] <- (si[, 1, ] > si[, 2, ]) %>% 
				tf$cast(tf$int32) %>% 
				tf$reduce_sum(0L)

			if (j %% 10 == 0){
				sprintf('test_accessibility | batch=%6.d/%6.d', j, length(bs)) %>%
					message()
			}
			
		}

		p <- S %>%
			tf$concat(axis = 0L) %>%
			tf$cast(tf$float32) %>% 
			tf$math$add(1L) %>% 
			tf$math$divide(sampling) %>% 
			tf$math$maximum(1/sampling) %>%
			tf$math$minimum(1 - 1/sampling) %>%
			as.numeric()

		if (all(contrasts[2:3][order(contrasts[2:3])] != contrasts[2:3]))
			p <- 1 - p

		bf <- log10(p) - log10(1 - p)

		res <- data.frame(bayes_factor_close = rep(NA, length(x)), bayes_factor_open = rep(NA, length(x)))
		res[include, ] <- data.frame(bayes_factor_close = bf, bayes_factor_open = -bf)

		mcols(x)[[sprintf('%s,%s', contrasts[2], contrasts[3])]] <- res
		x
	}
) # test_accessibility

