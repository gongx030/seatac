setOldClass('rpytools.call.RNNEncoderModel')
setOldClass('python.builtin.RNNEncoderModel')
setClassUnion('RNNEncoderModel', members = c('rpytools.call.RNNEncoderModel', 'python.builtin.RNNEncoderModel'))

RNNEncoderModel <- reticulate::PyClass(
	'RNNEncoderModel',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self, 
			kmer_embedding,
			kmers_size,
			gru_dim,
			bin_size,
			block_size,
			vae,
			rate = 0.1
		){
			super()$`__init__`()

			self$bin_size <- as.integer(bin_size)
			self$block_size <- as.integer(block_size)
			self$n_bins_per_block <- as.integer(self$block_size / self$bin_size)

			self$embedding <- tf$keras$layers$Embedding(as.integer(kmers_size), kmer_embedding)
			self$dropout_1 <- tf$keras$layers$Dropout(rate)

			self$bigru <- tf$keras$layers$Bidirectional(layer = tf$keras$layers$GRU(units = gru_dim))
			self$dropout_2 <- tf$keras$layers$Dropout(rate)
			self$dense <- tf$keras$layers$Dense(vae$encoder$latent_dim)

			self$vae <- vae
			self$vae$trainable <- FALSE

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

			x <- x %>% self$embedding()
			x <- x %>% self$dropout_1()
			x <- x %>% self$bigru()
			x <- x %>% self$dropout_2()
			z <- x %>% self$dense()
			y <- z %>% self$vae$decoder()
			list(
				latent = z,
				y = y$x_pred,
				nucleosome = y$nucleosome
			)
		}
	)
)




#' prepare_data
#'
setMethod(
	'prepare_data',
	signature(
		model = 'RNNEncoderModel',
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
			types = c('nucleoatac', 'full_nucleoatac')
		) 

		flog.info('prepare_data | get_latent_representation')
		d$latent_representation <- model$vae %>% get_latent_representation(d$vplots, batch_size = 128L)
		
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
		model = 'RNNEncoderModel',
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
		train_loss <- tf$keras$losses$MeanSquaredError(reduction = 'none')

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		# updating step
		train_step <- function(x, y){

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model(x)
				loss <- train_loss(y, res$latent) %>%
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
			loss <- train_loss(y, res$latent) %>%
			list(loss = loss)
		}

		test_step <- tf_function(test_step) # convert to graph mode

		for (epoch in seq_len(epochs)){

			loss_train <- NULL 
			iter <- make_iterator_one_shot(x$train)
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$kmers, batch$latent_representation)
				loss_train <- c(loss_train, as.numeric(res$loss))
			})

			loss_test <- NULL
			iter <- make_iterator_one_shot(x$test)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$kmers, batch$latent_representation)
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
		model = 'RNNEncoderModel',
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
		predicted_nucleosome <- matrix(0, length(x), x@n_bins_per_window)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting | batch=%4.d/%4.d', i, n_batch))

			b <- starts[i]:ends[i]

			d <- select_blocks(
				x[b],
				block_size = model$block_size,
				min_reads = 0,
				with_kmers = TRUE
			)

			res <- model(d$kmers)

			y_pred <- res$y %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, model$n_bins_per_block, 1L)) %>%
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


