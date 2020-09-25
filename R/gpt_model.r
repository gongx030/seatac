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
			num_layers, 
			d_model,
			num_heads, 
			n_intervals,
			bin_size,
			dff, 
			block_size,
			kmers_size,
			rate = 0.1
		){
			super()$`__init__`()

			self$d_model <- d_model
			self$num_layers <- num_layers
			self$n_intervals <- as.integer(n_intervals)
			self$bin_size <- as.integer(bin_size)
			self$block_size <- as.integer(block_size)
			self$n_bins_per_block <- as.integer(self$block_size / self$bin_size)
			self$embedding <- tf$keras$layers$Embedding(as.integer(kmers_size), d_model)

			self$pos_encoding <- positional_encoding(self$block_size, d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)

			self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(d_model, num_heads = num_heads, dff = dff, rate = rate))

			self$dropout <- tf$keras$layers$Dropout(rate)

			self$pool_1d <- tf$keras$layers$AveragePooling1D(self$bin_size, self$bin_size)

			self$final_layer <- tf$keras$layers$Dense(
				units = self$n_intervals, 
				activation = 'softmax'
			)

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

			x <- x %>% self$embedding()
			x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))
			x <- x + self$pos_encoding
			x <- self$dropout(x, training = training)
							    
			for (i in seq_len(self$num_layers)){
				x <- self$enc_layers[[i - 1]](x, training = training, mask = mask)	# zero-based
			}

			x <- self$pool_1d(x)
			x <- self$final_layer(x)
			x <- x %>% tf$transpose(c(0L, 2L, 1L))
			x <- x %>% tf$expand_dims(3L)
			x 
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

		d <- select_blocks(
			x,
			block_size = model$block_size,
			min_reads = min_reads,
			with_kmers = TRUE,
			with_predicted_counts = TRUE,
			types = c('nucleoatac', 'full_nucleoatac')
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
#		learning_rate <- CustomSchedule(model$d_model)
#		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		train_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')
		bce <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		# updating step
		train_step <- function(x, y){

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				y_pred <- model(x)
				loss <- train_loss(y, y_pred) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
					
			})

			gradients <- tape$gradient(loss, model$trainable_variables)
			list(gradients, model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()

			list(loss = loss)

		} # train_step

		train_step <- tf_function(train_step) # convert to graph mode

		test_step <- function(x, y){
			y_pred <- model(x)
			loss <- train_loss(y, y_pred) %>%
				tf$reduce_sum(shape(1L, 2L)) %>%
				tf$reduce_mean()
			list(loss = loss)
		}

		test_step <- tf_function(test_step) # convert to graph mode


		for (epoch in seq_len(epochs)){

			loss_train <- NULL 
			iter <- make_iterator_one_shot(x$train)
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$kmers, batch$predicted_vplots)
				loss_train <- c(loss_train, as.numeric(res$loss))
			})

			loss_test <- NULL
			iter <- make_iterator_one_shot(x$test)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- test_step(batch$kmers, batch$predicted_vplots)
				loss_test <- c(loss_test, as.numeric(res$loss))
			})

			flog.info(sprintf('epoch=%6.d/%6.d | train_loss=%13.7f | test_loss=%13.7f', epoch, epochs, mean(loss_train), mean(loss_test)))

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

		n_blocks_per_window <- as.integer(x@n_bins_per_window - model$n_bins_per_block + 1)
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

#			if (i == 1 || i %% 100 == 0)
#				flog.info(sprintf('predicting | batch=%4.d/%4.d', i, n_batch))

			b <- starts[i]:ends[i]

			xi <- x[b]$kmers %>%
				tf$cast(tf$int32) %>%
				tf$expand_dims(-1L) %>%
				tf$expand_dims(-1L) %>%
				tf$image$extract_patches(
					sizes = c(1L, model$block_size, 1L, 1L),
					strides = c(1L, x@bin_size, 1L, 1L),
					rates = c(1L, 1L, 1L, 1L),
					padding = 'VALID'
				) %>%
				tf$squeeze(axis = 2L)

			xi <- xi %>%
				tf$reshape(c(xi$shape[[1]] * xi$shape[[2]], model$block_size))

			yi <- model(xi)

			yi <- yi %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, model$n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			yi <- aperm(yi, c(1, 3, 2))
			dim(yi) <- c(length(b), x@n_bins_per_window * x@n_intervals)

			predicted_counts[b, ] <- yi

		}

		mcols(x)$predicted_counts <- predicted_counts

		x
	}
) # predict


