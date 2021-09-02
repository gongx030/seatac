#' VplotPredModel
#' 
#' A VAE model for V-plot
#'
#' @param latent_dim Latent dimension (default: 10L)
#' @param block_size Block size in base pairs (default: 640L)
#' @param bin_size Bin size in base pairs(default: 5L) 
#' @param fragment_size_range  Fragment size ranges (default: c(0L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
VplotPredModel <- function(
	decoder = NULL,
	latent_dim = 10L,
	bin_size = 5L,
	block_size = 640L,
	embedding_dim = 32L,
	kernel_size = 5L,
	num_blocks = 5L,
	fragment_size_range  = c(0L, 320L),
	fragment_size_interval = 10L,
	n_samples = 1L,
	ratio = 0.5,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$latent_dim <- latent_dim
		self$bin_size <- bin_size
		self$block_size <- block_size
		self$n_bins_per_block <- as.integer(block_size / bin_size)
		self$n_samples <- n_samples 
		self$fragment_size_range <- fragment_size_range
		self$fragment_size_interval <- fragment_size_interval
		br <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		self$n_intervals <- length(br) - 1L

		self$decoder <- decoder
		self$decoder$trainable <- FALSE

		self$sequence_embedding <- tf$keras$layers$Conv1D(filters = embedding_dim, kernel_size = bin_size, strides = bin_size, padding = 'valid')
		self$bin_embedding <- tf$keras$layers$Conv1D(filters = embedding_dim, kernel_size = kernel_size, strides = 1L, padding = 'same')
		self$positional_embedding <- tf$keras$layers$Embedding(input_dim = self$n_bins_per_block, output_dim = embedding_dim)

		self$fnet <- tf$keras$Sequential(lapply(1:num_blocks, function(i) FNetLayer(
			embedding_dim = embedding_dim,
			rate = rate
		)))

		self$conv <- tf$keras$layers$Conv1D(filters = self$n_intervals, kernel_size = 3L, strides = 1L, padding = 'same')

		self$dense <- tf$keras$layers$Dense(self$n_intervals)

		function(x, ..., training = TRUE){

			vplots_from <- x$vplots_from %>%
				tf$squeeze(3L) %>%
				tf$transpose(shape(0L, 2L, 1L)) 

			if (training){
				p <- (tf$random$uniform(shape(1L, self$n_bins_per_block, 1L)) < ratio) %>%
					tf$cast(tf$float32)
				vplots_from <- vplots_from * p + tf$zeros_like(vplots_from) * (1 - p)
			}

			vplots_from <- vplots_from %>% 
				self$bin_embedding()

			s <- x$sequence %>%
				tf$one_hot(5L) %>%
				self$sequence_embedding()
			
			pe <- tf$range(0L, self$n_bins_per_block) %>%
				self$positional_embedding()

			vplots_to <- (vplots_from + s + pe) %>%
				self$fnet() %>%
				self$conv() %>%
				tf$transpose(shape(0L, 2L, 1L)) %>%
				tf$expand_dims(3L) %>%
				tf$keras$activations$softmax(1L)

			list(
				vplots_to = vplots_to
			)
		}
	})
}


#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VplotPredModel object, initialized by `new('VplotPredModel', model = VplotPredModel(...))`
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
		model = 'VplotPredModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		samples = NULL,
		...
	){

		stopifnot(all(samples %in% x@samples))
		stopifnot(length(samples) >= 2L)

		stopifnot(!is.null(assays(x)[['predicted_counts']]))
		stopifnot(!is.null(rowData(x)[['vae_z_mean']]))

		vplots_from <- NULL
		vplots_to <- NULL
		sequence <- NULL

		for (i in 1:(length(samples) - 1L)){

			from <- rowData(x)$batch == samples[i]
			to <- rowData(x)$batch == samples[i + 1L]

			v_from <- assays(x[from])$predicted_counts %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>%
				tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L))

			v_to <- assays(x[to])$predicted_counts %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>%
				tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L))

			s <- rowData(x[from])$sequence %>%
				as.matrix() %>%
				factor(c('N', 'A', 'C', 'G', 'T')) %>%
				as.numeric() %>%
				matrix(nrow = sum(from), ncol = x@window_size) %>%
				tf$cast(tf$int64)
			s <- s - 1L

			vplots_from <- c(vplots_from, v_from)
			vplots_to <- c(vplots_to, v_to)
			sequence <- c(sequence, s)
		}

		vplots_from <- vplots_from %>% tf$concat(axis = 0L)
		vplots_to <- vplots_to %>% tf$concat(axis = 0L)
		sequence <- sequence %>% tf$concat(axis = 0L)

		list(vplots_from = vplots_from, vplots_to = vplots_to, sequence = sequence)
	}
)


#' fit
#'
#' Fit a VplotPredModel
#'
#' @param model a VplotPredModel object, initialized by `new('VplotPredModel', model = VplotPredModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @return a VplotPredModel
#'
setMethod(
	'fit',
	signature(
		model = 'VplotPredModel',
		x = 'tf_dataset'
	),
	function(
		 model,
		 x,
		 batch_size = 128L,
		 epochs = 100L,
		 learning_rate = 1e-3,
		 test_size = 0.15,
		 compile = TRUE
	 ){

		x <- x %>%
			split_dataset(test_size, batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')  # loss for the V-plot

		train_step <- function(batch){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(batch, training = TRUE)
				loss <- reconstrution_loss(batch$vplots_to, res$vplots_to)  %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
	 		})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		for (epoch in seq_len(epochs)){
			loss_train <- NULL
			iter <- x$train %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch)
				loss_train <- rbind(loss_train, sapply(res, as.numeric))
			})

			if (epoch == 1 || epoch %% 10 == 0){

				loss_train <- colMeans(loss_train)
				sprintf('epoch=%6.d/%6.d | train | %s', epoch, epochs, paste(sapply(1:length(loss_train), function(i) sprintf('%s=%13.7f', names(loss_train)[i], loss_train[i])), collapse = ' | ')) %>%
					message()

				loss_test <- NULL
				iter <- x$test %>%
					dataset_shuffle(1000L) %>%
					make_iterator_one_shot()
				res <- until_out_of_range({
					batch <- iterator_get_next(iter)
					res <- model@model(batch, training = FALSE)
					res <- list(
						loss = reconstrution_loss(batch$vplots_to, res$vplots_to)  %>%
							tf$reduce_sum(shape(1L, 2L)) %>%
							tf$reduce_mean()
					)
					loss_test <- rbind(loss_test, sapply(res, as.numeric))
				})

				loss_test <- colMeans(loss_test)
				sprintf('epoch=%6.d/%6.d |  test | %s', epoch, epochs, paste(sapply(1:length(loss_test), function(i) sprintf('%s=%13.7f', names(loss_test)[i], loss_test[i])), collapse = ' | ')) %>%
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
#' @param model a trained VplotPredModel object
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
		model = 'VplotPredModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L, # v-plot per batch
		samples = NULL,
		...
	){

		d <- model %>%
			prepare_data(x, samples = samples, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()

		from <- NULL
		to <- NULL
		to_pred <- NULL

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch, training = FALSE)
			from <- c(from, batch$vplots_from)
			to <- c(to, batch$vplots_to)
			to_pred <- c(to_pred, res$vplots_to)
		})

		to_pred <- to_pred %>% tf$concat(axis = 0L)


		modified_counts <- matrix(0, nrow = nrow(x), ncol = ncol(x))
		j <- unlist(lapply(1:(length(samples) - 1L), function(i) which(rowData(x)$batch == samples[i + 1L])))

		modified_counts[j, ] <- to_pred %>% 
			tf$reshape(shape(to_pred$shape[[1]], -1L)) %>%
			as.matrix()
		dimnames(modified_counts) <- dimnames(x)
		assays(x)$modified_counts <- modified_counts
		
		x
	}
)



#' predicted_fragment_size
#'
#' Predict fragment size at the center of the Vplot
#'
#' @param model a trained VplotPredModel object
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
		model = 'VplotPredModel',
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
