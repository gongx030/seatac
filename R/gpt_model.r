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

#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'TransformerDecoderModel',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		batch_size = 32L,
		batch_size_predict = 1L,
		batch_size_select = 128L,
		epochs = 100L,
		min_reads_per_block = 10,
		plot = FALSE
	){

		# Adopted from https://www.tensorflow.org/tutorials/text/transformer#optimizer
		learning_rate <- CustomSchedule(model$d_model)
		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		bce <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')

		# updating step
		train_step <- function(x, y){

			w <- tf$reduce_sum(y, 1L, keepdims = TRUE) > 0
			w <- w %>% tf$cast(tf$float32)
			w <- w %>% tf$squeeze(3L)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				y_pred <- model(x)

				loss <- (w * bce(y, y_pred))  %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()

			})

			gradients <- tape$gradient(loss, model$trainable_variables)
			list(gradients, model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()

			loss

		} # train_step

		train_step <- tf_function(train_step) # convert to graph mode

		inputs <- select_blocks(x, block_size = block_size, min_reads = min_reads_per_block, batch_size = batch_size_select, with_kmers = TRUE)

		dataset <- tf$data$Dataset$from_tensor_slices(inputs) %>%
			dataset_shuffle(1000) %>%
			dataset_repeat() %>%
			dataset_batch(batch_size)

		iter <- make_iterator_one_shot(dataset)

		for (epoch in seq_len(epochs)){

			batch <- iterator_get_next(iter)
			xi <- batch$kmers
			yi <- batch$vplot

			loss <- train_step(xi, yi)

			if (epoch %% 1000 == 0){

				# evaluating the predicted performance
				x_sample <- sample(x, 500)
				x_pred <- model %>% predict(x_sample, batch_size = batch_size_predict, with_kmers = TRUE)
				x_pred <- add_nucleosome_signal(x_pred)

				if (plot){
					vplot(x_pred, 'predicted_counts', main = sprintf('epoch=%d', epoch))
				}

				rmse_mnase_seatac <- sqrt(rowSums((x_pred$mnase_scaled - x_pred$nucleosome_signal)^2))
				rmse_mnase_nucleoatac <- sqrt(rowSums((x_pred$mnase_scaled - x_pred$nucleoatac_scaled)^2))

				if (!is.null(x_pred$full_nucleoatac_scaled)){
					rmse_nucleoatac_seatac <- sqrt(rowSums((x_pred$full_nucleoatac_scaled - x_pred$nucleosome_signal)^2))
					rmse_nucleoatac_nucleoatac <- sqrt(rowSums((x_pred$full_nucleoatac_scaled - x_pred$nucleoatac_scaled)^2))
				}else{
					rmse_nucleoatac_seatac <- 0
					rmse_nucleoatac_nucleoatac <- 0
				}

				flog.info(sprintf('evaluating | epoch=%6.d/%6.d | loss=%13.7f | rmse_mnase=%.3f(%.3f) | rmse_nucleoatac=%.3f(%.3f)',  epoch, epochs, loss, mean(rmse_mnase_seatac), mean(rmse_mnase_nucleoatac), mean(rmse_nucleoatac_seatac), mean(rmse_nucleoatac_nucleoatac)))

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


