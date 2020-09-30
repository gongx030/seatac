setOldClass('rpytools.call.CVaeModel')
setOldClass('python.builtin.CVaeModel')
setClassUnion('CVaeModel', members = c('rpytools.call.CVaeModel', 'python.builtin.CVaeModel'))

VplotEncoder <- reticulate::PyClass(
	'KmersEncoder',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			kmers_size
		){
			super()$`__init__`()

			self$kmers_size <- kmers_size

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

		}
	)
)


CVaeEncoder <- reticulate::PyClass(
	'CVaeEncoder',
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
CVaeDecoder <- reticulate::PyClass(
	'CVaeDecoder',
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
CVaeModel <- reticulate::PyClass(
	'CVaeModel',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			latent_dim,
			n_intervals,
			bin_size,
			block_size,
			filters0 = 64L,
			rate = 0.1,
			kmers_size
		){

			super()$`__init__`()

			if (block_size %% bin_size != 0)
				stop('block_size must be a multiple of bin_size')

			self$block_size <- block_size
			self$bin_size <- bin_size
			self$n_bins_per_block <- as.integer(block_size / bin_size)

			self$kmers_encoder <- KmersEncoder(
				kmers_size = kmers_size
			)

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
		call = function(self, x, g, training = TRUE){

			y <- g %>% self$kmers_encoder()
			browser()

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
		model = 'CVaeModel',
		x = 'VplotsKmers'
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
			with_kmers = TRUE,
			types = c('nucleoatac', 'mnase', 'full_nucleoatac')
		) %>%
			tensor_slices_dataset()
		d
	}
)




#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'CVaeModel',
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
		model = 'CVaeModel',
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
		model = 'CVaeModel',
		x = 'VplotsKmers'
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
