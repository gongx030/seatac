#' VaeModel
#' 
#' Build a VAE model for V-plot of multiple ATAC-seq datasets. This model takes the stacked V-plots of the same genomic regions as the input. 
#'
#' @param x a Vplots or a VplotsList object
#' @param latent_dim Latent dimension (default: 10L)
#' @param filters0 Filter size after the latent layer (default: 128L)
#' @param filters Initial filter size of convolution layer (default: 32L)
#' @param kernel_size Kernel size in convolution and deconvolution  layers (default: 3L)
#' @param downsample_layers Downsample layers (default: 4L)
#' @param upsample_layers Upsample layers (default: 4L)
#' @param strides Convolution strides 
#' @param momentum Momentum in BatchNormalization layer (default: 0.8)
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
#' @return a VaeModel object
#'
#' @export
#'
VaeModel <- function(
	x,
	latent_dim = 10L,
	filters0 = 128L,
	filters = 32L,
	kernel_size = 3L,
	downsample_layers = 4L,
	upsample_layers = 4L,
	strides = c(2L, 2L),
	momentum = 0.8,
	rate = 0.1,
	name = NULL
){

	if (!is(x, 'Vplots') && !is(x, 'VplotsList')){
		stop('x must be a Vplots or a VplotsList object')
	}

	if (!(is.integer(latent_dim) && length(latent_dim) == 1 & latent_dim >= 2)){
		stop('latent_dim must be an integer greater than or equal to 2')
	}

	if (!(is.integer(filters0) && length(filters0) == 1 & filters0 >= 1)){
		stop('filters0 must be a positive integer')
	}

	if (!(is.integer(filters) && length(filters) == 1 & filters >= 1)){
		stop('filters must be a positive integer')
	}

	if (!(is.integer(kernel_size) && length(kernel_size) == 1 & kernel_size >= 1)){
		stop('kernel_size must be a positive integer')
	}

	if (!(is.integer(downsample_layers) && length(downsample_layers) == 1 & downsample_layers >= 1)){
		stop('downsample_layers (layers of encoders) must be a positive integer')
	}

	if (!(is.integer(upsample_layers) && length(upsample_layers) == 1 & upsample_layers >= 1)){
		stop('upsample_layers (layers of decoders) must be a positive integer')
	}

	if (!(is.integer(strides) && length(strides) == 2 & upsample_layers >= 1)){
		stop('strides must be positive integers')
	}

	if (is(x, 'Vplots')){
		n_channels <- x@n_samples
	}else if (is(x, 'VplotsList')){
		n_channels <- 1L
	}

	model <- .build_VaeModel(
		n_samples = x@n_samples,
		n_channels = n_channels,
		latent_dim = latent_dim,
		block_size = x@window_size,
		bin_size = x@bin_size,
		filters0 = filters0,
		filters = filters,
		upsample_layers = upsample_layers,
		downsample_layers = downsample_layers,
		fragment_size_range = x@fragment_size_range,
		fragment_size_interval = x@fragment_size_interval
	)
	model <- new('VaeModel', model = model)
	model
}


.build_VaeModel <- function(
	n_samples = 1L,
	n_channels = 1L,
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
	strides = c(2L, 2L),
	momentum = 0.8,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$n_samples <- n_samples
		self$n_channels <- n_channels
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

		self$dense_1 <- tf$keras$layers$Dense(units = (self$n_channels + 1L) * self$latent_dim)

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
			loc  = tf$zeros(shape(self$n_channels, latent_dim)),
			scale_identity_multiplier = 1
		)

		function(x, ..., training = TRUE){

			batch_size <- x$batch$shape[[1]]

			# when x$vplots was originally sparse, it appears we need to specify the dimensions here
			# otherwise it will complain about unknown channel dimensions. 
			# this will only happen in compiled function (e.g. after the tf_function call)
			x$vplots <- x$vplots %>% 
				tf$reshape(shape(batch_size, self$n_intervals, self$n_bins_per_block, self$n_channels))

			fragment_size <- x$batch %>% 
				self$dense_fragment_size() %>%
				tf$transpose(shape(0L, 2L, 1L)) %>%
				tf$expand_dims(2L) 

			y <- (x$vplots + fragment_size) %>% 
				self$encoder() %>%
				self$dense_1() %>%
				tf$reshape(shape(batch_size, self$n_channels + 1L, self$latent_dim))

			q_m <- y[, 1:self$n_channels, , drop = FALSE]
			q_v <- tf$nn$softplus(y[, self$n_channels + 1L, , drop = FALSE] + 1e-3)

			posterior <- tfp$distributions$MultivariateNormalDiag(
				loc = q_m,
				scale_diag = q_v
			)

			if (training){
				z <- posterior$sample()
				b <- x$batch %>% tf$one_hot(self$n_channels)
			}else{
				z <- posterior$mean()
				b <- tf$zeros(shape(batch_size, self$n_channels), dtype = tf$int64) %>% tf$one_hot(self$n_channels)
			}

			x_pred <- list(z, b) %>% 
				tf$concat(2L) %>% 
				tf$reshape(shape(batch_size * self$n_channels, self$latent_dim + self$n_channels)) %>%
				self$decoder(training = training) %>%
				tf$reshape(shape(batch_size, self$n_channels, self$n_intervals, self$n_bins_per_block)) %>%
				tf$transpose(shape(0L, 2L, 3L, 1L)) %>%
				tf$keras$activations$softmax(1L)

			list(
				posterior = posterior, 
				z = z, 
				vplots = x_pred
			)
		}
	})
}


#' validate
#'
#' Validate the input data for current model
#'
#' @param model a VaeModel object
#' @param x a Vplots object
#'
#' @return a logical value of whether the data are valid.
#'
#' @export
#'
setMethod(
	'validate',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x
	){
		stopifnot(model@model$n_samples == x@n_samples)
		stopifnot(!is.null(assays(x, withDimnames = FALSE)$counts))
		stopifnot(all(model@model$fragment_size_range == x@fragment_size_range))
		stopifnot(all(model@model$fragment_size_interval == x@fragment_size_interval))
		stopifnot(all(model@model$bin_size == x@bin_size))
		stopifnot(all(model@model$block_size == x@window_size))
		return(TRUE)

	}
)


#' validate
#'
#' Validate the input data for current model
#'
#' @param model a VaeModel object
#' @param x a VplotsList object
#'
#' @return a logical value of whether the data are valid.
#'
#' @export
#'
setMethod(
	'validate',
	signature(
		model = 'VaeModel',
		x = 'VplotsList'
	),
	function(
		model,
		x
	){
		stopifnot(model@model$n_samples == x@n_samples)
		for (i in seq_len(x@n_samples)){
			stopifnot(!is.null(assays(x[[i]], withDimnames = FALSE)$counts))
			stopifnot(x[[i]]@n_samples == 1L)
		}
		stopifnot(all(model@model$fragment_size_range == x@fragment_size_range))
		stopifnot(all(model@model$fragment_size_interval == x@fragment_size_interval))
		stopifnot(all(model@model$bin_size == x@bin_size))
		stopifnot(all(model@model$block_size == x@window_size))
		return(TRUE)

	}
)


#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VaeModel object, initialized by `VaeModel()`
#' @param x a Vplots object
#' @param ... Other arguments
#'
#' @return a list that include `vplots` and `batch`
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
		...
	){

		validate(model, x)

		v <- assays(x, withDimnames = FALSE)$counts %>%
			summary()

		v <- tf$sparse$SparseTensor(
			indices = v[, 1:2] %>% as.matrix() %>% tf$cast(tf$int64) - 1L,
			values = v[, 3] %>% tf$cast(tf$float32),
			dense_shape = shape(dim(x)['grange'],  dim(x)['sample'] * dim(x)['interval'] * dim(x)['bin'])
		) %>%
			tf$sparse$reshape(dim(x)) %>%
			tf$sparse$transpose(shape(0L, 2L, 3L, 1L)) %>%
			tf$sparse$reorder() 

		batch <- tf$range(dim(x)['sample']) %>%
			tf$reshape(shape(1L, dim(x)['sample'])) %>%
			tf$`repeat`(nrow(x), axis = 0L)

		list(vplots = v, batch = batch)
	}
)


#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VaeModel object, initialized by `VaeModel()`
#' @param x a VplotsList object
#' @param ... Other arguments
#'
#' @return a list that include `vplots` and `batch`
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
		x,
		...
	){

		validate(model, x)

		vplots <- NULL
		batch <- NULL
		offset <- 0L

		for (i in seq_len(x@n_samples)){

			v <- assays(x[[i]], withDimnames = FALSE)$counts %>%
				summary()

			v <- tf$sparse$SparseTensor(
				indices = v[, 1:2] %>% as.matrix() %>% tf$cast(tf$int64) - 1L,
				values = v[, 3] %>% tf$cast(tf$float32),
				dense_shape = shape(dim(x[[i]])['grange'],  dim(x[[i]])['sample'] * dim(x[[i]])['interval'] * dim(x[[i]])['bin'])
			) %>%
				tf$sparse$reshape(dim(x[[i]])) %>%
				tf$sparse$transpose(shape(0L, 2L, 3L, 1L)) %>%
				tf$sparse$reorder() 

			b <- tf$range(x[[i]]@n_samples) %>%
				tf$reshape(shape(1L, x[[i]]@n_samples)) %>%
				tf$`repeat`(nrow(x[[i]]), axis = 0L)

			b <- b + offset
			offset <- offset + x[[i]]@n_samples
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
#' @param beta Beta sequences (default: 1)
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
		 beta = 1
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
			res <- until_out_of_range2({
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


#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a list
#' @param ... Additional arguments to fit('VaeModel`, 'tf_dataset', ...)
#'
#' @export
#' @return a VaeModel
#'
setMethod(
	'fit',
	signature(
		model = 'VaeModel',
		x = 'list'
	),
	function(
		 model,
		 x,
		 ...
	 ){
		x <- x %>% tensor_slices_dataset()
		fit(model, x, ...)
	}
)


#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a Vplots object
#' @param ... Additional arguments to fit('VaeModel`, 'list', ...)
#'
#' @export
#' @return a VaeModel
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
		 ...
	 ){
		x <- model %>% prepare_data(x)
		fit(model, x, ...)
	}
)

#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a VplotsList object
#' @param ... Additional arguments to fit('VaeModel`, 'list', ...)
#'
#' @export
#' @return a VaeModel
#'
setMethod(
	'fit',
	signature(
		model = 'VaeModel',
		x = 'VplotsList'
	),
	function(
		 model,
		 x,
		 ...
	 ){
		x <- model %>% prepare_data(x)
		fit(model, x, ...)
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
#' @param nucleosome Whether or not return predicted nucleosome as rowData(x)$predicted_nucleosome
#' @param fragment_size_threshold Fragment size threshold for nucleosome reads (default: 150L)
#' @param ... Additional arguments
#' @importFrom SummarizedExperiment assays<-
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
		vplots = FALSE,
		nucleosome = TRUE,
		fragment_size_threshold = 150L,
		...
	){

		validate(model, x)

		res <- .predict_counts(
			model, 
			x,
			batch_size = batch_size, # v-plot per batch
			vplots = vplots,
			nucleosome = nucleosome ,
			fragment_size_threshold = fragment_size_threshold
		)

		dimnames(res$z)[1:2] <- list(rownames(x), x@dimdata$sample$name)
		dimnames(res$z_stddev)[1:2] <- list(rownames(x), x@dimdata$sample$name)
		rowData(x)[['vae_z_mean']] <- res$z
		rowData(x)[['vae_z_stddev']] <- res$z_stddev

		if (vplots){
			# dimnames(res$predicted_vplots[j, , drop = FALSE]) <- dimnames(x[[i]])
			assays(x, withDimnames = FALSE)$predicted_counts <- res$predicted_vplots
		}

		if (nucleosome){
			rowData(x)[['predicted_nucleosome']] <- res$predicted_nucleosome
		}
		x
	}
)

#' predict
#'
#' Predict counts
#'
#' @param model a trained VaeModel object
#' @param x a VplotsList object
#' @param batch_size Batch size (default: 256L)
#' @param vplots Whether or not return predicted Vplots as assays(x)$predicted_counts
#' @param nucleosome Whether or not return predicted nucleosome as rowData(x)$predicted_nucleosome
#' @param fragment_size_threshold Fragment size threshold for nucleosome reads (default: 150L)
#' @param ... Additional arguments
#' @importFrom SummarizedExperiment assays<-
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
		x = 'VplotsList'
	),
	function(
		model,
		x,
		batch_size = 256L, # v-plot per batch
		vplots = FALSE,
		nucleosome = TRUE,
		fragment_size_threshold = 150L,
		...
	){

		validate(model, x)

		res <- .predict_counts(
			model, 
			x,
			batch_size = batch_size, # v-plot per batch
			vplots = vplots,
			nucleosome = nucleosome ,
			fragment_size_threshold = fragment_size_threshold
		)

		start <- 1L
		for (i in seq_len(length(x))){
			j <- start:(start + length(x[[i]]) - 1)
			zi <- res$z[j, , , drop = FALSE]
			zi_stddev <- res$z_stddev[j, , , drop = FALSE]
			dimnames(zi)[1:2] <- list(rownames(x[[i]]), x[[i]]@dimdata$sample$name)
			dimnames(zi_stddev)[1:2] <- list(rownames(x[[i]]), x[[i]]@dimdata$sample$name)
			rowData(x[[i]])[['vae_z_mean']] <- zi
			rowData(x[[i]])[['vae_z_stddev']] <- zi_stddev

			if (vplots){
				# dimnames(res$predicted_vplots[j, , drop = FALSE]) <- dimnames(x[[i]])
				assays(x[[i]], withDimnames = FALSE)$predicted_counts <- res$predicted_vplots[j, , drop = FALSE]
			}

			if (nucleosome){
				rowData(x[[i]])[['predicted_nucleosome']] <- res$predicted_nucleosome[j, , , drop = FALSE]
			}

			start <- start + length(x[[i]])

		}
		x
	}
)

.predict_counts <- function(
	model, 
	x, 
	batch_size = 256L, # v-plot per batch
	vplots = FALSE,
	nucleosome = TRUE,
	fragment_size_threshold = 150L,
	...
){
	stopifnot(is.numeric(fragment_size_threshold) && fragment_size_threshold >= x@fragment_size_range[1] && fragment_size_threshold <= x@fragment_size_range[2])

	iter <- model %>% prepare_data(x) %>%
		tensor_slices_dataset() %>%
		dataset_batch(batch_size) %>%
		make_iterator_one_shot()

	z <- NULL
	z_stddev <- NULL

	predicted_vplots <- NULL
	predicted_nucleosome  <- NULL

	if (nucleosome){
		if (is(x, 'Vplots')){
			is_nucleosome <- x@dimdata$interval$center >= fragment_size_threshold
		}else if (is(x, 'VplotsList')){
			is_nucleosome <- x[[1]]@dimdata$interval$center >= fragment_size_threshold
		}
	}

	res <- until_out_of_range2({
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
		if (nucleosome){
			predicted_nucleosome <- c(
				predicted_nucleosome, 
				res$vplots %>%
					tf$boolean_mask(is_nucleosome, 1L) %>%
					tf$reduce_sum(1L)
			)
		}
	})

	z <- z %>% tf$concat(axis = 0L)
	z_stddev <- z_stddev %>% tf$concat(axis = 0L)

	z <- z %>%  as.array()
	z_stddev  <- z_stddev  %>% as.array()

	if (vplots){

		predicted_vplots <- predicted_vplots %>% 
			tf$concat(axis = 0L) %>%
			tf$transpose(shape(0L, 3L, 1L, 2L)) 

		predicted_vplots <- predicted_vplots %>% 
			tf$reshape(shape(dim(predicted_vplots)[1], prod(dim(predicted_vplots)[-1]))) %>%
			as.matrix()

	}

	if (nucleosome){
		predicted_nucleosome <- predicted_nucleosome %>%
			tf$concat(axis = 0L) %>%
			tf$transpose(shape(0L, 2L, 1L))  %>%
			as.array()
	}
	list(
		z = z,
		z_stddev = z_stddev,
		predicted_vplots = predicted_vplots,
		predicted_nucleosome = predicted_nucleosome
	)
}


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

		stopifnot(is.numeric(width) && width >= 0 && width <= x@window_size)
		validate(model, x)

		d <- model %>%
			prepare_data(x, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()

		is_center <- x@dimdata$bin$position >= -width / 2 & x@dimdata$bin$position <= width / 2

		fragment_size <- NULL
		predicted_fragment_size <- NULL

		res <- until_out_of_range2({
			batch <- iterator_get_next(iter)
			batch$vplots <- batch$vplots %>% 
				tf$sparse$to_dense() %>%
				scale01()
			res <- model@model(batch, training = FALSE)
			fs <- batch$vplots %>% tf$boolean_mask(is_center, 2L) %>% tf$reduce_sum(2L) 
			w <- fs %>% tf$reduce_sum(1L, keepdims = TRUE)
			fs <- fs / tf$where(w > 0, w, tf$ones_like(w))
			fs_pred <- res$vplots %>% tf$boolean_mask(is_center, 2L) %>% tf$reduce_mean(2L) 
			fragment_size <- c(fragment_size, fs)
			predicted_fragment_size <- c(predicted_fragment_size, fs_pred)
		})

		fragment_size <- fragment_size %>% 
			tf$concat(axis = 0L) %>%
			tf$transpose(shape(0L, 2L, 1L)) %>%
			as.array()
		predicted_fragment_size <- predicted_fragment_size %>% 
			tf$concat(axis = 0L) %>%
			tf$transpose(shape(0L, 2L, 1L)) %>%
			as.array()
		rowData(x)[['fragment_size']] <- fragment_size
		rowData(x)[['predicted_fragment_size']] <- predicted_fragment_size
		x
	}
)


#' load_model
#'
#' Load a pretrained VaeModel
#'
#' @param model a VaeModel object
#' @param dir Model directory. The index file should be 'filename.index' and the data file should be 'filename.data-00000-of-00001'
#'
#' @return a VaeModel object
#'
#' @export
#'
setMethod(
	'load_model',
	signature(
		model = 'VaeModel'
	),
	function(
		model,
		dir = 'character'
	){


		model_index_file <- sprintf('%s.index', dir)
		model_data_file <- sprintf('%s.data-00000-of-00001', dir)

		stopifnot(file.exists(model_index_file))
		stopifnot(file.exists(model_data_file))

		vplots <- tf$random$uniform(shape(1L, model@model$n_intervals, model@model$n_bins_per_block, model@model$n_samples))
		batch <- tf$zeros(shape(1L, model@model$n_samples), dtype = tf$int64)
		res <- model@model(list(vplots = vplots, batch = batch))
		model@model$load_weights(dir)
		model
	}
)

