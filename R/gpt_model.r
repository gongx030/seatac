setOldClass('rpytools.call.TransformerDecoderModel')
setOldClass('python.builtin.TransformerDecoderModel')
setClassUnion('TransformerDecoderModel', members = c('rpytools.call.TransformerDecoderModel', 'python.builtin.TransformerDecoderModel'))

CustomSchedule <- reticulate::PyClass(
	'CustomSchedule',
	inherit = tf$keras$optimizers$schedules$LearningRateSchedule,
	list(
		`__init__` = function(
			self, 
			d_model,
			warmup_steps = 4000
		){
			super()$`__init__`()

			self$d_model <- tf$cast(d_model, tf$float32)
			self$warmup_steps <- warmup_steps

			NULL
		},
		`__call__` = function(self, step){

			step <- tf$cast(step, tf$float32)
			arg1 <- tf$math$rsqrt(step)
			arg2 <- step * (self$warmup_steps ^ (-1.5))
		    
			tf$math$rsqrt(self$d_model) * tf$math$minimum(arg1, arg2)
		}
	)
)


TransformerDecoderModel <- reticulate::PyClass(
	'TransformerDecoderModel',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self, 
			d_model,
			num_layers, 
			num_heads, 
			dff, 
			kmers_size,
			vae,
			rate = 0.1,
			step_size = 2L
		){
			super()$`__init__`()

			self$d_model <- d_model

			stopifnot(self$d_model %% num_heads == 0)

			self$num_layers <- num_layers

			self$vae <- vae
			self$vae$trainable <- FALSE

			self$pos_encoding <- positional_encoding(self$vae$block_size, self$d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)

			self$mask <- tf$constant(1:self$vae$block_size %in% seq(1, self$vae$block_size, by = step_size), tf$bool)

			self$embedding <- tf$keras$layers$Embedding(as.integer(kmers_size), self$d_model)
			self$dropout_1 <- tf$keras$layers$Dropout(rate)

			self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(self$d_model, num_heads = num_heads, dff = dff, rate = rate))

			self$dropout_2 <- tf$keras$layers$Dropout(rate)
			self$dense <- tf$keras$layers$Dense(units = model_vae$encoder$latent_dim)

			NULL

		},
		call = function(self, x, training = TRUE){

			x <- x %>% self$embedding()
			x <- x %>% tf$boolean_mask(self$mask, axis = 1L)
			x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))
			x <- x + tf$boolean_mask(self$pos_encoding, self$mask, axis = 1L)
			x <- x %>% self$dropout_1()
							    
			for (i in seq_len(self$num_layers)){
				x <- self$enc_layers[[i - 1]](x)	# zero-based
			}

			x <- x %>% tf$reduce_mean(1L)
			x <- x %>% self$dropout_2()
			z <- x %>% self$dense()
			res  <- z %>% self$vae$decoder()
			res
		}
	)
)




#' prepare_data
#'
setMethod(
	'prepare_data',
	signature(
		model = 'TransformerDecoderModel',
		x = 'VplotsKmers'
	),
	function(
		 model,
		 x,
		 min_reads = 10
	 ){

		flog.info(sprintf('prepare_data | min_reads=%d', min_reads))

		d <- select_blocks(
			x,
			block_size = model$vae$block_size,
			min_reads = min_reads,
			with_kmers = TRUE,
			types = c('nucleoatac', 'full_nucleoatac')
		) 

		flog.info(sprintf('prepare_data | number of samples=%d', d$vplots$shape[[1]]))

		d <- d %>%
			tensor_slices_dataset()
		d
	}
)



#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'TransformerDecoderModel',
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

		# Adopted from https://www.tensorflow.org/tutorials/text/transformer#optimizer
		learning_rate <- CustomSchedule(model$d_model)
		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		train_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		# updating step
		train_step <- function(x, y, w){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model(x)
				loss <- (tf$squeeze(w, 3L) * train_loss(res$x_pred, y)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
					
			})
			gradients <- tape$gradient(loss, model$trainable_variables)
			list(gradients, model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(loss = loss)
		} # train_step

		test_step <- function(x, y, w){
			res <- model(x)
			loss <- (tf$squeeze(w, 3L) * train_loss(res$x_pred, y)) %>%
				tf$reduce_sum(shape(1L, 2L)) %>%
				tf$reduce_mean()
			list(loss = loss)
		}

		train_step <- tf_function(train_step) # convert to graph mode
		test_step <- tf_function(test_step) # convert to graph mode

		for (epoch in seq_len(epochs)){

			loss_train <- NULL 
			iter <- make_iterator_one_shot(x$train)
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$kmers, batch$vplots, batch$weight)
				loss_train <- c(loss_train, as.numeric(res$loss))
			})

			loss_test <- NULL
			iter <- make_iterator_one_shot(x$test)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- test_step(batch$kmers, batch$vplots, batch$weight)
				loss_test <- c(loss_test, as.numeric(res$loss))
			})

			flog.info(sprintf('epoch=%6.d/%6.d | train_loss=%13.7f | test_loss=%13.7f', epoch, epochs, mean(loss_train), mean(loss_test)))

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
	'predict',
	signature(
		model = 'TransformerDecoderModel',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		batch_size = 1L, # v-plot per batch
		with_kmers = FALSE
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		n_blocks_per_window <- as.integer(x@n_bins_per_window - model$vae$n_bins_per_block + 1)
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)
		predicted_nucleosome <- matrix(0, length(x), x@n_bins_per_window)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting | batch=%4.d/%4.d', i, n_batch))

			b <- starts[i]:ends[i]

			d <- select_blocks(
				x[b],
				block_size = model$vae$block_size,
				min_reads = 0,
				with_kmers = TRUE
			)

			res <- model(d$kmers)

			y_pred <- res$x_pred %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, model$vae$n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			y_pred <- aperm(y_pred, c(1, 3, 2))
			dim(y_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)
			predicted_counts[b, ] <- y_pred

			predicted_nucleosome[b, ] <- res$nucleosome %>%
				tf$expand_dims(2L) %>%
				tf$expand_dims(3L) %>%
				tf$transpose(shape(0L, 2L, 1L, 3L)) %>%
				tf$reshape(c(length(b), n_blocks_per_window, 1L, model$vae$n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks()  %>%
				tf$squeeze(shape(1L, 3L)) %>%
				as.matrix()

		}

		mcols(x)$predicted_counts <- predicted_counts
		mcols(x)$predicted_nucleosome <- predicted_nucleosome

		x
	}
) # predict

