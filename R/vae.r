#' VaeModel
#' 
#' A VAE model for V-plot
#'
#' @param latent_dim Latent dimension (default: 10L)
#' @param block_size Block size in base pairs (default: 640L)
#' @param bin_size Bin size in base pairs(default: 5L) 
#' @param filters0 Filter size after the latent layer (default: 128L)
#' @param filters Initial filter size of convolution layer (default: 32L)
#' @param kernel_size Kernel size in convolution and deconvolution  layers (default: 3L)
#' @param downsample_layers Downsample layers (default: 4L)
#' @param upsample_layers Upsample layers (default: 4L)
#' @param fragment_size_range  Fragment size ranges (default: c(0L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param n_batches Number of batches (default: 1L)
#' @param strides Convolution strides 
#' @param momentum Momentum in BatchNormalization layer (default: 0.8)
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
	 filters = 32L,
	 kernel_size = 3L,
	 downsample_layers = 4L,
	 upsample_layers = 4L,
	 fragment_size_range  = c(0L, 320L),
	 fragment_size_interval = 10L,
	 n_batches = 1L,
	 strides = c(2L, 2L),
	 momentum = 0.8,
	 rate = 0.1,
	 name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$latent_dim <- latent_dim
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
		self$n_batches <- n_batches

		stopifnot(self$n_intervals %% strides[1]^upsample_layers == 0)
		stopifnot(self$n_bins_per_block %% strides[2]^upsample_layers == 0)

		self$conv <- lapply(1:downsample_layers, function(i) tf$keras$Sequential(list(
			tf$keras$layers$Conv2D(
				filters = filters,
				kernel_size = kernel_size,
				strides = strides
			),
			tf$keras$layers$BatchNormalization(momentum = momentum),
			tf$keras$layers$Activation('relu')
		))) %>%
			tf$keras$Sequential()

		self$encoder <- tf$keras$Sequential(list(
			self$conv,
			tf$keras$layers$Flatten()
		))

		self$dense_1 <- tf$keras$layers$Dense(units = 2 * self$latent_dim)

		self$dense_fragment_size <- tf$keras$layers$Embedding(n_batches,  self$n_intervals)

		interval_dim0 <- as.integer(self$n_intervals / strides[1]^upsample_layers)
		window_dim0 <- as.integer(self$n_bins_per_block / strides[2]^upsample_layers)
		output_dim0 <- as.integer(window_dim0 * interval_dim0 * filters0)

		self$deconv <- lapply(1:upsample_layers, function(i) tf$keras$Sequential(list(
			tf$keras$layers$Conv2DTranspose(
				filters = filters, 
				kernel_size = kernel_size,
				strides = strides,
				padding = 'same'
			),
			tf$keras$layers$Activation('relu'),
			tf$keras$layers$BatchNormalization(momentum = momentum)
		))) %>%
			tf$keras$Sequential()

		self$decoder <- tf$keras$Sequential(list(
			tf$keras$layers$Dense(units = output_dim0,activation = 'relu'),
			tf$keras$layers$Reshape(target_shape = c(interval_dim0, window_dim0, filters0)),
			self$deconv,
			tf$keras$layers$Conv2D(
				filters = 1L,
				kernel_size = 1L,
				padding = 'same'
			)
		))

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(latent_dim)),
			scale_identity_multiplier = 1
		)

		function(x, ..., training = TRUE){

			batch_size <- x$batch$shape[[1]]

			# when x$vplots was originally sparse, it appears we need to specify the dimensions here
			# otherwise it will complain about unknown channel dimensions. 
			# this will only happen in compiled function (e.g. after the tf_function call)
			x$vplots <- x$vplots %>% 
				tf$reshape(shape(batch_size, self$n_intervals, self$n_bins_per_block, 1L))

			fragment_size <- x$batch %>% 
				self$dense_fragment_size() %>%
				tf$expand_dims(2L) %>%
				tf$expand_dims(3L)

			y <- (x$vplots + fragment_size) %>% 
				self$encoder() %>%
				self$dense_1()

			q_m <- y[, 1:self$latent_dim]
			q_v <- tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3)

			posterior <- tfp$distributions$MultivariateNormalDiag(
				loc = q_m,
				scale_diag = q_v
			)

			if (training){
				z <- posterior$sample()
				b <- x$batch %>% tf$one_hot(self$n_batches)
			}else{
				z <- posterior$mean()
				b <- tf$zeros(shape(z$shape[[1]]), dtype = tf$int64) %>% tf$one_hot(self$n_batches)
			}

			x_pred <- list(z, b) %>% 
				tf$concat(1L) %>% 
				self$decoder(training = training) %>%
				tf$keras$activations$softmax(1L)

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
#' @param ... Other arguments
#'
#' @return a list that include `vplots` and `batch`
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
		x,
		...
	){

		d <- list()

		vplots <- NULL
		batch <- NULL

		for (i in 1:x@n_samples){

			j <- colData(x)$batch == x@samples[i]
			v <- assays(x)$counts[, j, drop = FALSE] %>%
				summary()

			v <- tf$sparse$SparseTensor(
				 indices = v[, 1:2] %>% as.matrix() %>% tf$cast(tf$int64) - 1L,
				 values = v[, 3] %>% tf$cast(tf$float32),
				 dense_shape = shape(nrow(x), x@n_intervals * x@n_bins_per_window)
			) %>%
				tf$sparse$reshape(shape(nrow(x), x@n_intervals, x@n_bins_per_window, 1L)) %>%
				tf$sparse$reorder()

			b <- rep(x@samples[i], nrow(x)) %>%
				factor(x@samples) %>%
				as.numeric() %>%
				tf$cast(tf$int64) 
			b <- b - 1L

			vplots <- c(vplots, v)
			batch <- c(batch, b)
		}

		vplots <- tf$sparse$concat(0L, vplots)
		batch <- batch %>% tf$concat(axis = 0L)

		list(vplots = vplots, batch = batch)
	}
)


#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 128L)
#' @param epochs Number of training epochs (default: 100L)
#' @param learning_rate Learning rate (default: 1e-3)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#' @param beta Beta sequences (default: 5e-5)
#'
#' @export
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
		 batch_size = 128L,
		 epochs = 100L,
		 learning_rate = 1e-3,
		 compile = TRUE,
		 beta = 5e-5
	 ){

		if (length(beta) == 1)
			beta <- rep(beta, epochs)
		beta <- tf$cast(beta, tf$float32)

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		train_step <- function(batch, b){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				batch$vplots <- batch$vplots %>% 
					tf$sparse$to_dense() %>%
					scale01()
				res <- model@model(batch, training = TRUE)
				loss_reconstruction <- reconstrution_loss(batch$vplots, res$vplots) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss_kl <- b * loss_kl
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

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		for (epoch in seq_len(epochs)){
			loss <- NULL
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch, beta[epoch])
				loss <- rbind(loss, sapply(res, as.numeric))
			})

			loss <- colMeans(loss)
			sprintf('epoch=%6.d/%6.d | beta=%.3e | %s', epoch, epochs, beta[epoch], paste(sapply(1:length(loss), function(i) sprintf('%s=%13.7f', names(loss)[i], loss[i])), collapse = ' | ')) %>%
				message()
		}
		model
	}
)

#' predict
#'
#' Predict counts
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#' @param vplots Whether or not return predicted Vplots as assays(x)$predicted_counts
#' @param ... Additional arguments
#'
#' @return a Vplots object 
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
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
		batch_size = 256L, # v-plot per batch
		vplots = TRUE,
		...
	){


		iter <- model %>%
			prepare_data(x, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size) %>%
			make_iterator_one_shot()

		z <- NULL
		z_stddev <- NULL

		if (vplots){
			predicted_vplots <- NULL
		}

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			batch$vplots <- batch$vplots %>% 
				tf$sparse$to_dense() %>%
				scale01()
			res <- model@model(batch, training = FALSE)
			z <- c(z, res$posterior$mean())
			z_stddev <- c(z_stddev, res$posterior$stddev())
			if (vplots){
				predicted_vplots <- c(predicted_vplots, res$vplots)
			}
		})

		z <- z %>% tf$concat(axis = 0L)
		z_stddev <- z_stddev %>% tf$concat(axis = 0L)

		z <- z %>% 
			tf$reshape(shape(x@n_samples, nrow(x), model@model$latent_dim)) %>%
			tf$transpose(shape(1L, 0L, 2L))

		z_stddev  <- z_stddev  %>% 
			tf$reshape(shape(x@n_samples, nrow(x), model@model$latent_dim)) %>%
			tf$transpose(shape(1L, 0L, 2L))

		rowData(x)[['vae_z_mean']] <- as.array(z)
		rowData(x)[['vae_z_stddev']] <- as.array(z_stddev)

		if (vplots){

			predicted_vplots <- predicted_vplots %>% 
				tf$concat(axis = 0L) %>%
				tf$reshape(shape(x@n_samples, nrow(x), x@n_intervals * x@n_bins_per_window)) %>%
				tf$transpose(shape(1L, 0L, 2L)) %>%
				tf$reshape(shape(nrow(x), x@n_samples * x@n_intervals * x@n_bins_per_window)) %>%
				as.matrix()

			dimnames(predicted_vplots) <- dimnames(x)
			assays(x)$predicted_counts <- predicted_vplots
		}

		x
	}
)



#' predicted_fragment_size
#'
#' Predict fragment size at the center of the Vplot
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#' @param width Central width (in bp) for computing fragment size disstribution (default: 100L)
#' @param ... Other arguments passed to prepare_data. 
#'
#' @return a Vplots object 
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'predict_fragment_size',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L, # v-plot per batch
		width = 100L,
		...
	){

		d <- model %>%
			prepare_data(x, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()

		is_center <- x@positions >= -width / 2 & x@positions <= width /2

		fragment_size <- NULL
		predicted_fragment_size <- NULL

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			batch$vplots <- batch$vplots %>% 
				tf$sparse$to_dense() %>%
				scale01()
			res <- model@model(batch, training = FALSE)
			fs <- batch$vplots %>% tf$boolean_mask(is_center, 2L) %>% tf$reduce_sum(2L) %>% tf$squeeze(2L)
			w <- fs %>% tf$reduce_sum(1L, keepdims = TRUE)
			fs <- fs / tf$where(w > 0, w, tf$ones_like(w))
			fs_pred <- res$vplots %>% tf$boolean_mask(is_center, 2L) %>% tf$reduce_mean(2L) %>% tf$squeeze(2L)
			fragment_size <- c(fragment_size, fs)
			predicted_fragment_size <- c(predicted_fragment_size, fs_pred)
		})

		fragment_size <- fragment_size %>% tf$concat(axis = 0L)
		predicted_fragment_size <- predicted_fragment_size %>% tf$concat(axis = 0L)
		rowData(x)[['fragment_size']] <- as.matrix(fragment_size)
		rowData(x)[['predicted_fragment_size']] <- as.matrix(predicted_fragment_size)
		x
	}
)

