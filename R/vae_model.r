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

		function(x, ...){

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
	
		function(x, training = TRUE, ...){

			y <- x %>%
				self$vplot_decoder(training = training)

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

		function(x, fragment_size, training = TRUE){

			posterior <- x %>% self$encoder()

			if (training){
				z <- posterior$sample()
			}else{
				z <- posterior$mean()
			}
			c <- list(z, fragment_size) %>% tf$concat(axis = 1L)
			x_pred <- c %>% self$decoder(training = training)

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
			matrix(length(x), 1L) %>% 
			tf$cast(tf$int32) %>%
			tf$math$subtract(1L) %>%
			tf$one_hot(x@n_samples) %>% 
			tf$squeeze(1L)

		fragment_size <- assays(x)$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window)) %>%
			tf$reduce_sum(2L) # fragment size per sample

		fragment_size <- batch %>% 
			tf$matmul(fragment_size, transpose_a = TRUE) %>%
			scale01()	# average fragment size per batch

		fragment_size <- batch %>%
			tf$matmul(fragment_size)	# proporate the batch-wise fragment size to each sample

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
			fragment_size = fragment_size
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

		sprintf('# of training samples: %d', length(x)) %>% 
			message()

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

		if (warmup > epochs)
			warmup <- epochs

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		train_step <- function(x, w, fragment_size, beta){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, fragment_size, training = TRUE)
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
				res <- train_step(batch$vplots, batch$weight, batch$fragment_size, beta[epoch])
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

#' predict_vplots
#' 
#' Predict V-plots
#' 
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#' @return a new Vplots object which assays have the predicted counts
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#' 
setMethod(
	'predict_vplots',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L # v-plot per batch
	){

		stopifnot(model@model$block_size == x@window_size)
		is_nucleosome <- model@model$centers >= 180 & model@model$centers <= 247

		d <- model %>% 
			prepare_data(x) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()
		x_pred <- NULL
		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch$vplots, batch$fragment_size, training = FALSE)
			x_pred <- c(x_pred, res$vplots)
		})
		x_pred <- x_pred %>% tf$concat(axis = 0L)

		assays(x)$predicted_counts <- x_pred %>% 
			tf$reshape(shape(length(x), -1L)) %>%
			as.matrix()
		x
	}
)


#' predict_nucleosome
#' 
#' Predict nucleosome signal
#' 
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#' @return a new Vplots object which assays have the predicted counts
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#' 
setMethod(
	'predict_nucleosome',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L # v-plot per batch
	){

		stopifnot(model@model$block_size == x@window_size)
		is_nucleosome <- model@model$centers >= 180 & model@model$centers <= 247

		d <- model %>% 
			prepare_data(x) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()
		nucleosome_signal <- NULL
		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch$vplots, batch$fragment_size, training = FALSE)

			ns_pred <- res$vplots %>% 
				tf$boolean_mask(is_nucleosome, axis = 1L) %>%
				tf$reduce_sum(1L, keepdims = TRUE)  %>%
				tf$squeeze(shape(1L, 3L)) %>%
				get_nucleosome_score()

			nucleosome_signal <- c(nucleosome_signal , ns_pred)
		})
		nucleosome_signal <- nucleosome_signal %>% tf$concat(axis = 0L)
		rowData(x)$nucleosome_signal <- nucleosome_signal %>% as.matrix()
		x
	}
)



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
		n0 <- length(include)	# number of total unique intervals
		n <- sum(include)	# number of qualitified genomic region to test
		include <- matrix(include, nrow = n0, ncol = x@n_samples, dimnames = list(NULL, x@samples))
		include[, !colnames(include) %in% contrasts[2:3]] <- FALSE
		include <- c(include)

		# there might be a memory issue if length(x) is too large
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
				model@model$decoder(training = FALSE) %>% 
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

		if (which(x@samples == contrasts[2]) > which(x@samples == contrasts[3]))
			p <- 1 - p

		bf <- log10(p) - log10(1 - p)

		res <- data.frame(bayes_factor_close = rep(NA, length(x)), bayes_factor_open = rep(NA, length(x)))
		res[include, ] <- data.frame(bayes_factor_close = bf, bayes_factor_open = -bf)

		mcols(x)[[sprintf('%s,%s', contrasts[2], contrasts[3])]] <- res
		x
	}
) # test_accessibility


#' test_nucleosome
#'
#' Test the nucleosome at the center of the block
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param sampling Number of sampling used to computing Bayes factor (default: 100L)
#' @param batch_size Batch size (default: 4L)
#'
#' @return a Vplots object that includes the nucleosome prediction (at the center of the block) in rowData fields (rowData(x)$nucleosome_score)
#' 				* mu: mean nucleosome scores from resampling
#' 				* std: std of nucleosome scores from resampling
#' 				* nu: sum of two shape parameters of beta distribution (nu = alpha + beta)
#' 				* alpha, beta: estimated shape parameter
#' 				* prob: probability of nucleosome score greater than 0.5, estimated by a beta distribution
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
#' 
setMethod(
	'test_nucleosome',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		sampling = 100L,
		batch_size = 4L
	){

    stopifnot(model@model$block_size == x@window_size)

		is_center <- block_center(model@model$block_size / model@model$bin_size)
    is_nucleosome <- model@model$centers >= 180 & model@model$centers <= 247

		d <- prepare_data(model, x) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()

		nucleosome <- NULL

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			n <- batch$vplots$shape[[1]]
			posterior <- model@model$encoder(batch$vplots)
			z <- posterior$sample(sampling)
			fragment_size <- batch$fragment_size %>%
				tf$expand_dims(0L) %>%
				tf$tile(shape(sampling, 1L, 1L))

			x_center <- list(z, fragment_size) %>%
				tf$concat(axis = 2L) %>%
				tf$reshape(shape(sampling * n, -1L)) %>%
				model@model$decoder(training = FALSE) %>%
				tf$boolean_mask(is_center, axis = 2L)

			ns <- x_center %>% 
				tf$boolean_mask(is_nucleosome, axis = 1L) %>%
				tf$reduce_sum(1L, keepdims = TRUE)  %>%
				tf$reduce_mean(2L) %>%
				tf$squeeze(shape(1L, 2L)) %>%
				get_nucleosome_score() %>%
				tf$reshape(shape(sampling, n)) 

			nucleosome <- c(nucleosome, ns)

		})

		nucleosome <-  nucleosome %>% tf$concat(axis = 1L)
		nucleosome <- nucleosome %>% 
			tf$transpose(shape(1L, 0L)) %>%
			as.matrix()

		rowData(x)$nucleosome_score <- data.frame(
			mu = rowMeans(nucleosome),
		  std = rowSds(nucleosome)
		) %>%
			mutate(nu = mu * (1 - mu) / std^2 - 1) %>%
		 	mutate(alpha = mu * nu, beta = (1 - mu) * nu) %>%
		 	mutate(prob = pbeta(0.5, shape1 = alpha, shape2 = beta, lower.tail = FALSE))
		x
	}
)

