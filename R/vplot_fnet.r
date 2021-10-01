#' VplotFnetModel
#' 
#' A FNet model for predicting the V-plot
#'
#' @param embedding_dim Embedding dimension (default: 32L)
#' @param kernel_size Kernel_size for input Vplot
#' @param num_blocks Number of FNet blocks
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
VplotFnetModel <- function(
	vae = NULL,
	embedding_dim = 32L,
	kernel_size = 5L,
	num_blocks = 5L,
	rate = 0.1,
	kmer = 5L,
	vplots_input = TRUE,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$vae <- vae
		self$vae$trainable <- FALSE

		self$kmer <- kmer
		self$vplots_input <- vplots_input

		if (vplots_input){
			self$bin_embedding <- tf$keras$layers$Dense(embedding_dim)
		}

		self$positional_embedding <- tf$keras$layers$Embedding(input_dim = vae$n_bins_per_block, output_dim = embedding_dim)

		self$sequence_embedding <- tf$keras$Sequential(list(
			tf$keras$layers$Embedding(as.integer(4^vae$bin_size), embedding_dim),
			tf$keras$layers$MaxPooling1D(vae$bin_size)
		))

		self$fnet <- tf$keras$Sequential(lapply(1:num_blocks, function(i) FNetLayer(
			embedding_dim = embedding_dim,
			rate = rate
		)))

		self$dense_fs <- tf$keras$layers$Dense(1L, activation = 'tanh')

		function(x, ..., training = TRUE){

			pe <- tf$range(0L, self$vae$n_bins_per_block) %>%
				self$positional_embedding()

			s <- x$sequence %>%
				self$sequence_embedding()

			if (self$vplots_input){

				b <- tf$zeros(shape(x$z$shape[[1]]), dtype = tf$int64) %>% tf$one_hot(self$vae$n_batches)

				if (training){
					prob <- tfp$distributions$MultivariateNormalDiag(
						loc = x$z,
						scale_diag = x$z_stddev
					)
					z <- prob$sample()
				}else{
					z <- x$z
				}

				v <- list(z, b) %>%
					tf$concat(1L) %>%
					self$vae$decoder(training = FALSE) %>%
					tf$keras$activations$softmax(1L) %>%
					tf$squeeze(3L) %>%
					tf$transpose(shape(0L, 2L, 1L)) %>%
					self$bin_embedding()

				x <- v + pe + s
			}else{
				x <- pe + s
			}

			nucleosome_diff  <- x %>%
				self$fnet() %>%
				self$dense_fs() %>%
				tf$squeeze(2L)

			list(
				nucleosome_diff = nucleosome_diff
			)
		}
	})
}


#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VplotFnetModel object, initialized by `new('VplotFnetModel', model = VplotFnetModel(...))`
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
		model = 'VplotFnetModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		samples = NULL,
		fragment_size_threshold = 150L,
		training = TRUE,
		...
	){

		stopifnot(all(samples %in% x@samples))
		stopifnot(length(samples) >= 2L)

		stopifnot(!is.null(rowData(x)[['sequence']]))
		stopifnot(!is.null(rowData(x)[['vae_z_mean']]))

		sequences <- NULL
		nucleosome_diff <- NULL
		z <- NULL
		z_stddev <- NULL

		kmers <- Reduce('paste0', expand.grid(lapply(1:x@bin_size, function(j)  c('A', 'C', 'G', 'T'))))
		nucleotides <- c('A', 'C', 'G', 'T')

		for (i in 1:(length(samples) - 1L)){

			from <- rowData(x)$batch == samples[i]
			to <- rowData(x)$batch == samples[i + 1L]

			s <- rowData(x[from])$sequence %>%
				as.matrix()
			prob <- table(factor(c(s), nucleotides))
			invalid <- !s %in% nucleotides
			if (any(invalid)){
				s[invalid] <- sample(nucleotides, sum(invalid), replace = TRUE, prob = prob / sum(prob))
			}

			# add randomly padding to the end of the sequences
			s <- cbind(s, sample(nucleotides, sum(from) * (model@model$kmer - 1), replace = TRUE) %>% matrix(sum(from), model@model$kmer - 1))

			s <- do.call('cbind', lapply(1:x@window_size, function(j){
				Reduce('paste0', as.data.frame(s[, j:(j + model@model$kmer - 1)])) %>%
					factor(kmers) %>%
					as.numeric()
			}))
			s <- s %>% tf$cast(tf$int64)
			s <- s - 1L

			if (training){
				nucleosome_diff <- c(nucleosome_diff, get_nucleosome_diff(x[from], x[to], fragment_size_threshold))
				z_stddev <- c(z_stddev, rowData(x[from])$vae_z_stddev %>% tf$cast(tf$float32))
			}

			z <- c(z, rowData(x[from])$vae_z_mean %>% tf$cast(tf$float32))
			sequences <- c(sequences, s)
		}

		sequences <- sequences %>% tf$concat(axis = 0L)
		z <- z %>% tf$concat(axis = 0L)

		
		if (training){
			z_stddev <- z_stddev %>% tf$concat(axis = 0L)
			nucleosome_diff <- nucleosome_diff %>% tf$concat(axis = 0L)
			list(
				z = z,
				z_stddev = z_stddev,
				sequence = sequences,
				nucleosome_diff = nucleosome_diff
			)
		}else{
			list(
				z = z,
				sequence = sequences
			)
		}

	}
)


#' fit
#'
#' Fit a VplotFnetModel
#'
#' @param model a VplotFnetModel object, initialized by `new('VplotFnetModel', model = VplotFnetModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @return a VplotFnetModel
#'
setMethod(
	'fit',
	signature(
		model = 'VplotFnetModel',
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

		reconstrution_loss <- tf$keras$losses$MeanSquaredError()

		train_step <- function(batch){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(batch, training = TRUE)
				loss <- reconstrution_loss(batch$nucleosome_diff, res$nucleosome_diff)

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

			if (epoch == 1 || epoch %% 1 == 0){

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
					loss <- reconstrution_loss(batch$nucleosome_diff, res$nucleosome_diff)
					res <- list(
						loss = loss
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
#' @param model a trained VplotFnetModel object
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
		model = 'VplotFnetModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L, # v-plot per batch
		contrast = NULL,
		observed = TRUE,
		...
	){

		field <- contrast[1]
		control <- contrast[2]
		treatment <- contrast[3]
		x <- x[rowData(x)[[field]] %in% c(treatment, control)]

		d <- model %>%
			prepare_data(x, samples = c(control, treatment), training = observed, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)


		iter <- d %>% make_iterator_one_shot()

		nucleosome_diff <- NULL
		nucleosome_diff_pred <- NULL

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch, training = FALSE)
			nucleosome_diff_pred <- c(nucleosome_diff_pred, res$nucleosome_diff)
			if (observed){
				nucleosome_diff <- c(nucleosome_diff, batch$nucleosome_diff)
			}
		})

		nucleosome_diff_pred <- nucleosome_diff_pred %>% tf$concat(axis = 0L)

		res <- granges(x[rowData(x)[[field]] == control])
	  mcols(res) <- mcols(res)
		mcols(res)$nucleosome_diff_pred <- as.matrix(nucleosome_diff_pred)

		if (observed){
			nucleosome_diff <- nucleosome_diff %>% tf$concat(axis = 0L)
			mcols(res)$nucleosome_diff <- as.matrix(nucleosome_diff)
		}
		res
	}
)

