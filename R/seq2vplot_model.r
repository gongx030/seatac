#' Seq2VplotModel
#' 
#' @export
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'Seq2VplotModel',
	signature(
		model = 'VaeModel'
	),
	function(
		model,
		d_model,
		num_layers, 
		num_heads, 
		dff, 
		kmers_size,
		step_size = 2L,
		rate = 0.1,
		name = NULL
	){

		keras_model_custom(name = name, function(self){

			self$d_model <- d_model

			stopifnot(self$d_model %% num_heads == 0)

			self$num_layers <- num_layers
			self$vae <- model@model	# a kerastools.model.RModel object
			self$vae$trainable <- FALSE

			self$pos_encoding <- positional_encoding(self$vae$block_size, self$d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)

			self$annotation_encoding_1 <- tf$keras$layers$Dense(8L, activation = 'relu')
			self$annotation_encoding_2 <- tf$keras$layers$Dense(self$d_model)

			self$embedding <- tf$keras$layers$Embedding(as.integer(kmers_size), self$d_model)
			self$dropout_1 <- tf$keras$layers$Dropout(rate)
		
			if (num_layers > 0){
				self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(self$d_model, num_heads = num_heads, dff = dff, rate = rate))
			}
	
			self$dense_1 <- tf$keras$layers$Dense(units = self$vae$decoder$vplot_decoder$vplot_height, activation = 'softmax')
			self$pool <- tf$keras$layers$MaxPooling1D(self$vae$bin_size, self$vae$bin_size)
	
			function(x, h, training = TRUE){

				x <- x %>% self$embedding()
				x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))

				annotation <- h %>% 
					self$annotation_encoding_1() %>%
					self$annotation_encoding_2()

				x <- x + self$pos_encoding + annotation
				x <- x %>% self$dropout_1()
		
				if (num_layers > 0){
					for (i in seq_len(self$num_layers)){
						x <- self$enc_layers[[i - 1]](x)	# zero-based
					}
				}

				x <- x %>% self$dense_1()
				x <- x %>% self$pool()
				x <- x %>% tf$transpose(shape(0L, 2L, 1L))
				x <- x %>% tf$expand_dims(3L)
				posterior <- self$vae$encoder(x)
				z <- posterior$mean()
				res <- self$vae$decoder(z)
				res$z <- z
				res
			}	
		})
	}
)

#' prepare_data
#' 
#' Prepare data for training/testing Seq2VplotModel model
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'prepare_data',
	signature(
		model = 'Seq2VplotModel',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		annotation
	){

		d <- list()

		y <- assays(x)$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32)

		w <- y %>% tf$reduce_sum(shape(1L), keepdims = TRUE)	# sum of reads per bin
		y <- y / tf$where(w > 0, w, tf$ones_like(w))	# scale so that the sum of each bar is one (softmax)

		d$z <- new('VaeModel', model = model@model$vae) %>% 
			encode(y, batch_size = 256L)	# for blocks

		d$kmers <- rowData(x)$kmers %>%
			tf$cast(tf$int32)

		g <- lapply(annotation, function(a){ 
			rowData(x)[[a]] %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(2L)
		})
		g <- tf$concat(g, axis = 2L)
		d$annotation <- g

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
		model = 'Seq2VplotModel',
		x = 'tf_dataset'
	),
	function(
		model,
		x,
		batch_size = 32L,
		epochs = 100L,
		test_size = 0.15,
		learning_rate = 1e-4
	){

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		train_loss <- tf$keras$losses$MeanSquaredError(reduction = 'none')

		x <- x %>% 
			dataset_shuffle(1000L) %>%
			split_dataset(test_size = test_size, batch_size = batch_size)

		# updating step
		train_step <- function(x, y, h){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				 res <- model@model(x, h)
				 loss <- train_loss(res$z, y) %>%
				 tf$reduce_mean()
			})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(loss = loss, predicted_vplots = res$vplots)
		} # train_step

		test_step <- function(x, y, h){
			res <- model@model(x, h)
			loss <- train_loss(res$z, y) %>%
			tf$reduce_mean()
			list(loss = loss, predicted_vplots = res$vplots)
		}

		train_step <- tf_function(train_step) # convert to graph mode
		test_step <- tf_function(test_step) # convert to graph mode

		for (epoch in seq_len(epochs)){
			loss_train <- NULL 
			iter <- make_iterator_one_shot(x$train)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$kmers, batch$z, batch$annotation)
				loss_train <- c(loss_train, as.numeric(res$loss))
			})
			loss_test <- NULL
			iter <- make_iterator_one_shot(x$test)
			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- test_step(batch$kmers, batch$z, batch$annotation)
				loss_test <- c(loss_test, as.numeric(res$loss))
			})
			message(sprintf('%s | fit | epoch=%6.d/%6.d | train_loss=%13.7f | test_loss=%13.7f', Sys.time(), epoch, epochs, mean(loss_train), mean(loss_test)))
		}
		model
	}
)

