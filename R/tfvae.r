#' tfVaeModel
#' 
#' A transcriptional factor VAE model for analyzing multiple ATAC-seq datasets
#'
#' @param latent_dim Latent dimension (default: 10L)
#' @param block_size Block size in base pairs (default: 640L)
#' @param bin_size Bin size in base pairs(default: 5L) 
#' @param filters0 Filter size after the latent layer (default: 128L)
#' @param fragment_size_range  Fragment size ranges (default: c(0L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
tfVaeModel <- function(
	 latent_dim = 10L,
	 block_size = 640L,
	 bin_size = 5L,
	 filters0 = 128L,
	 fragment_size_range  = c(0L, 320L),
	 fragment_size_interval = 10L,
	 rate = 0.1,
	 name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

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
		
		self$encoder <- VplotEncoder(
			latent_dim = latent_dim,
			filters = c(32L, 32L, 32L),
			kernel_size = c(3L, 3L, 3L),
			window_strides = c(2L, 2L, 2L),
			interval_strides = c(2L, 2L, 2L),
			rate = rate
		)

		self$decoder <- VplotDecoder(
			vplot_width = self$n_bins_per_block,
			vplot_height = self$n_intervals,
			filters0 = filters0,
			filters = c(32L, 32L, 1L),
			kernel_size = c(3L, 3L, 3L),
			window_strides = c(2L, 2L, 2L),
			interval_strides = c(2L, 2L, 1L),
			rate = rate
		)

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(latent_dim)),
			scale_identity_multiplier = 1
		)

		function(x, training = TRUE){

			browser()
			posterior <- x %>% self$encoder()

			if (training){
				z <- posterior$sample()
			}else{
				z <- posterior$mean()
			}
			c <- list(z, fragment_size) %>% tf$concat(axis = 1L)
			x_pred <- c %>% self$decoder(training = training)

			list(
				posterior = posterior, 
				z = z, 
				vplots = x_pred
			)
		}
	})
}


#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a tfVaeModel object, initialized by `new('tfVaeModel', model = tfVaeModel(...))`
#' @param x a Vplots object
#'
#' @return a list that include `vplots`, `weight` and `batch`
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'prepare_data',
	signature(
		model = 'tfVaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_field = 'batch',
		motif_field = 'motif',
		...
	){

		d <- list()

		stopifnot(!is.null(rowData(x)[[batch_field]]))
		stopifnot(!is.null(rowData(x)[[motif_field]]))

		batch <- rowData(x)$batch %>%
			factor(x@samples) %>%
			as.numeric() %>%
			tf$cast(tf$int32) %>%
			tf$math$subtract(1L) %>%
			tf$one_hot(x@n_samples) 

		vplots <- assays(x)$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
			scale_vplot()

		w <- tf$reduce_sum(vplots, 1L, keepdims = TRUE) > 0
		w <- w %>% tf$squeeze(3L)
		w <- w %>% tf$cast(tf$float32)

		motif <- rowData(x)[[motif_field]] %>%
			as.matrix() %>%
			tf$cast(tf$float32)

		list(vplots = vplots, batch = batch, weight = w, motif = motif)
	}
)


#' fit
#'
#' Fit a tfVaeModel
#'
#' @param model a tfVaeModel object, initialized by `new('tfVaeModel', model = tfVaeModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @return a tfVaeModel
#'
setMethod(
	'fit',
	signature(
		model = 'tfVaeModel',
		x = 'tf_dataset'
	),
	function(
		 model,
		 x,
		 batch_size = 32L,
		 epochs = 100L,
		 learning_rate = 1e-4,
		 warmup = 50L,
		 compile = TRUE
	 ){

		if (warmup > epochs)
			warmup <- epochs

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		train_step <- function(x, w, fragment_size, beta){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, fragment_size, training = TRUE)
				loss_reconstruction <- (w * reconstrution_loss(x, res$vplots)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction + beta * loss_kl
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

		beta <- c(seq(0, 1, length.out = warmup), rep(1, epochs - warmup))

		for (epoch in seq_len(epochs)){
			loss <- NULL 
			loss_reconstruction <- NULL
			loss_kl <- NULL
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch$vplots, batch$weight, batch$fragment_size, beta[epoch])
				loss <- c(loss, as.numeric(res$loss))
				loss_reconstruction <- c(loss_reconstruction, as.numeric(res$loss_reconstruction))
				loss_kl <- c(loss_kl, as.numeric(res$loss_kl))
			})

			sprintf('epoch=%6.d/%6.d | beta=%9.3f | recon_loss=%13.7f | kl_loss=%13.7f | loss=%13.7f', epoch, epochs, beta[epoch], mean(loss_reconstruction), mean(loss_kl), mean(loss)) %>%
				message()

		}
		model
	}
)
