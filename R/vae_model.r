setOldClass('rpytools.call.VaeModel')
setOldClass('python.builtin.VaeModel')

#' @export
setClassUnion('VaeModel', members = c('rpytools.call.VaeModel', 'python.builtin.VaeModel'))

VplotEncoder <- PyClass(
	'VplotEncoder',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			filters = c(32L, 32L, 32L),
			kernel_size = c(3L, 3L, 3L),
			window_strides = c(2L, 2L, 2L),
			interval_strides = c(2L, 2L, 2L)
		){
			super()$`__init__`()

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

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

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
	)
)

VplotDecoder <- PyClass(
	'VplotDecoder',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self, 
			vplot_width,	# width of the image
			vplot_height, # height of the image
			filters0 = 64L,
			filters = c(32L, 32L, 1L),
			kernel_size = c(3L, 3L, 3L),
			window_strides = c(2L, 2L, 2L),
			interval_strides = c(2L, 2L, 1L)
		){
			super()$`__init__`()

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

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

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
	)
)


#' Decoder


VaeEncoder <- PyClass(
	'VaeEncoder',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self, 
			latent_dim,
			rate = 0.1
		){
			super()$`__init__`()

			self$latent_dim <- latent_dim

			self$vplot_encoder <- VplotEncoder()

			self$dense_1 <- tf$keras$layers$Dense(
				units = 2 * self$latent_dim
			)

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

			y <- x %>%
				self$vplot_encoder() %>%
				self$dense_1()

			tfp$distributions$MultivariateNormalDiag(
				loc = y[, 1:self$latent_dim],
				scale_diag = tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3),
				name = 'posterior'
			)

		}
	)
)


#' Decoder
#' 
#' A transformer decoder to recover an image
#'
VaeDecoder <- PyClass(
	'VaeDecoder',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self, 
			vplot_width,	# width of the image
			vplot_height, # height of the image
			filters0 = 64,
			rate = 0.1
		){

			super()$`__init__`()

			self$vplot_decoder <- VplotDecoder(
				filters0 = filters0,																								 
				vplot_width = vplot_width, 
				vplot_height = vplot_height
			)

			self$final_layer <- tf$keras$layers$Dense(
				units = 1L,
				activation = 'sigmoid'
			)

			NULL
		},
		call = function(self, x, training = TRUE){

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
	)
)


#' VaeModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
VaeModel <- PyClass(
	'VaeModel',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			latent_dim = 10L,
			n_intervals,
			bin_size = 5L,
			block_size,
			filters0 = 128L,
			rate = 0.1
		){

			super()$`__init__`()

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


			NULL
		},
		call = function(self, x, training = TRUE){

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
	)
)



#' prepare_data
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
		min_reads = 10
	){

		d <- select_blocks(
			x,
			block_size = model$block_size,
			min_reads = min_reads,
			with_kmers = FALSE,
			types = c('nucleoatac', 'mnase', 'full_nucleoatac')
		) %>%
#			dataset_cache() %>% # https://www.tensorflow.org/guide/data_performance#reducing_memory_footprint
			tensor_slices_dataset()
		d
	}
)




#' fit
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
		checkpoint_dir = NULL
	){

		optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		train_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')
		bce <- tf$keras$losses$BinaryCrossentropy()

		if (!is.null(checkpoint_dir)){
			ckpt <- tf$train$Checkpoint(model = model, optimizer = optimizer)
			ckpt_manager <- tf$train$CheckpointManager(ckpt, checkpoint_dir, max_to_keep = 5)
		}

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		train_step <- function(x, w, y){

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				res <- model(x)

				loss_reconstruction <- (tf$squeeze(w, 3L) * train_loss(x, res$x_pred)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()

				loss_kl <- (res$posterior$log_prob(res$z) - model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()

				loss_nucleosome <- train_loss(y, res$nucleosome) %>%
					tf$reduce_mean()

				loss <- loss_reconstruction + loss_kl + 100 * loss_nucleosome

			})

			gradients <- tape$gradient(loss, model$trainable_variables)
			list(gradients, model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()

			list(
				loss = loss,
				loss_reconstruction = loss_reconstruction,
				loss_kl = loss_kl,
				loss_nucleosome = loss_nucleosome
			)
		}
		train_step <- tf_function(train_step) # convert to graph mode

		test_step <- function(x, w, y){

			res <- model(x)
			metric_nucleoatac <- bce(y, res$nucleosome)
			metric_test <- (tf$squeeze(w, 3L) * train_loss(x, res$x_pred)) %>%
				tf$reduce_sum(shape(1L, 2L)) %>%
				tf$reduce_mean()
			list(
				metric_test = metric_test,
				metric_nucleoatac = metric_nucleoatac
			)
		}
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

			flog.info(sprintf('epoch=%6.d/%6.d | train_recon_loss=%13.7f | train_kl_loss=%13.7f | train_nucleosome_loss=%13.7f | train_loss=%13.7f | test_recon_loss=%13.7f | test_nucleosome=%13.7f', epoch, epochs, mean(loss_train_reconstruction), mean(loss_train_kl), mean(loss_train_nucleosome), mean(loss_train), mean(metric_test), mean(metric_nucleoatac)))

			if (epoch %% 5 == 0 && !is.null(checkpoint_dir)){
				ckpt_save_path <- ckpt_manager$save()
				flog.info(sprintf('saving checkpoint at %s', ckpt_save_path))
			}
		}

		model
	}
)

#'
setMethod(
	'evaluate',
	signature(
		model = 'VaeModel',
		x = 'tf_dataset'
	),
	function(
		model,
		x,
		batch_size = 256L
	){

		bce <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')

		predict_step <- function(x, y, y_full){
			posterior <- model$encoder(x)
			res <- posterior$mean() %>% model$decoder()
			metric_nucleosome <- bce(y_full, res$nucleosome)
			metric_nucleoatac <- bce(y_full, y)
			mnase <- tf$reduce_mean(batch$mnase, 1L)
			list(
				mnase = mnase,
				metric_nucleosome = metric_nucleosome,
				metric_nucleoatac = metric_nucleoatac
			)
		}
		predict_step <- tf_function(predict_step) # convert to graph mode

		x <- x %>% 
			dataset_batch(batch_size)

		df <- NULL
		iter <- make_iterator_one_shot(x)
		until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- predict_step(batch$vplots, batch$nucleoatac, batch$full_nucleoatac)
			df <- rbind(df, data.frame(
				n = as.numeric(batch$n),
				mnase = as.numeric(res$mnase),
				metric_nucleosome = as.numeric(res$metric_nucleosome),
				metric_nucleoatac = as.numeric(res$metric_nucleoatac)
			))
		})
		df
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

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		n_blocks_per_window <- as.integer(x@n_bins_per_window - model$n_bins_per_block + 1)
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)
		predicted_nucleosome <- matrix(0, length(x), x@n_bins_per_window)

		predict_step <- function(x){
			posterior <- x %>% model$encoder()
			res <- posterior$mean() %>% model$decoder()
			res
		}
		predict_step <- tf_function(predict_step) # convert to graph mode

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting | batch=%4.d/%4.d', i, n_batch))

			b <- starts[i]:ends[i]

			d <- select_blocks(
				x[b],
				block_size = model$block_size,
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
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, model$n_bins_per_block, 1L)) %>%
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
				tf$reshape(c(length(b), n_blocks_per_window, 1L, model$n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks()  %>%
				tf$squeeze(shape(1L, 3L)) %>%
				as.matrix()
		}

		mcols(x)$predicted_counts <- predicted_counts
		mcols(x)$predicted_nucleosome <- predicted_nucleosome

		x
	}
) # predict
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

		starts <- seq(1, x$shape[[1]], by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > x$shape[[1]]] <- x$shape[[1]]
		n_batch <- length(starts)

		res <- list()

		for (i in 1:n_batch){
			b <- starts[i]:ends[i]
			res[[i]] <- model(x[b, , , ])
		}

		x_pred <- tf$concat(lapply(res, function(r) r$x_pred), axis = 0L)
		nucleosome <- tf$concat(lapply(res, function(r) r$nucleosome), axis = 0L)
		z <- tf$concat(lapply(res, function(r) r$z), axis = 0L)

		list(x_pred = x_pred, nucleosome = nucleosome, z = z)
	}
)

