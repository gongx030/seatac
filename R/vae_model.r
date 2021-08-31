#' VaeModel
#' 
#' A VAE model for V-plot
#'
#' @param latent_dim Latent dimension (default: 10L)
#' @param block_size Block size in base pairs (default: 640L)
#' @param bin_size Bin size in base pairs(default: 5L) 
#' @param filters0 Filter size after the latent layer (default: 128L)
#' @param fragment_size_range  Fragment size ranges (default: c(0L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
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
	 n_samples = 1L,
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
		self$n_samples <- n_samples

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

		self$dense_fragment_size <- tf$keras$layers$Embedding(n_samples,  self$n_intervals)


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
				b <- x$batch %>% tf$one_hot(self$n_samples)
			}else{
				z <- posterior$mean()
				b <- tf$zeros(shape(z$shape[[1]]), dtype = tf$int64) %>% tf$one_hot(self$n_samples)
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
#' @param weight Whether or not include positional weight
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
		x,
		...
	){

		d <- list()

		vplots <- assays(x)$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
			scale_vplot()

		batch <- rowData(x)$batch %>%
			factor(x@samples) %>%
			as.numeric() %>%
			tf$cast(tf$int64) 
		batch <- batch - 1L

		list(vplots = vplots, batch = batch)
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

		train_step <- function(x, b){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, training = TRUE)
				loss_reconstruction <- reconstrution_loss(x$vplots, res$vplots) %>%
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

			if (epoch == 1 || epoch %% 10 == 0){
				loss <- colMeans(loss)
				sprintf('epoch=%6.d/%6.d | beta=%.3e | %s', epoch, epochs, beta[epoch], paste(sapply(1:length(loss), function(i) sprintf('%s=%13.7f', names(loss)[i], loss[i])), collapse = ' | ')) %>%
					message()
			}
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
#'
#' @return a Vplots object with the mean and stddev of the latent representation
#'  (reducedDim(x, 'vae_z_mean') and reducedDim(x, 'vae_z_stddev'))
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
		reduction = 'vae_z_mean',
		reduction_stddev = 'vae_z_stddev',
		vplots = TRUE,
		...
	){

		d <- model %>%
			prepare_data(x, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()
		latent <- NULL
		latent_stddev <- NULL

		if (vplots){
			predicted_vplots <- NULL
		}

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch, training = FALSE)
			latent <- c(latent, res$z)
			latent_stddev <- c(latent_stddev, res$posterior$stddev())
			if (vplots){
				predicted_vplots <- c(predicted_vplots, res$vplots)
			}
		})

		latent <- latent %>% tf$concat(axis = 0L)
		latent_stddev <- latent_stddev %>% tf$concat(axis = 0L)
		rowData(x)[[reduction]] <- as.matrix(latent)
		rowData(x)[[reduction_stddev]] <- as.matrix(latent_stddev)

		if (vplots){
			predicted_vplots <- predicted_vplots %>% 
				tf$concat(axis = 0L) %>%
				tf$reshape(shape(nrow(x), -1L)) %>% 
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
#'
#' @return a Vplots object with the mean and stddev of the latent representation
#'  (reducedDim(x, 'vae_z_mean') and reducedDim(x, 'vae_z_stddev'))
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
