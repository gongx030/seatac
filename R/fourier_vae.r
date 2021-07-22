#' FourierVaeModel
#' 
#' A VAE model for V-plot, with fNet like encoder
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
FourierVaeModel <- function(
	latent_dim = 10L,
	block_size = 640L,
	bin_size = 5L,
	filters0 = 128L,
	filters = 32L,
	kernel_size = 3L,
	downsample_layers = 4L,
	upsample_layers = 4L,
	fragment_size_range  = c(0L, 320L),
	fragment_size_interval = 10L,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$latent_dim <- latent_dim
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
			downsample_layers = downsample_layers,
			filters = filters,
			kernel_size = kernel_size,
			rate = rate
		)

		self$encoder_f <- VplotEncoder(
			downsample_layers = downsample_layers,
			filters = filters,
			kernel_size = kernel_size,
			rate = rate
		)

		self$dense_1 <- tf$keras$layers$Dense(units = 2 * self$latent_dim)

		self$decoder <- VplotDecoder(
			vplot_width = self$n_bins_per_block,
			vplot_height = self$n_intervals,
			filters0 = filters0,
			upsample_layers = upsample_layers,
			filters = filters,
			kernel_size = kernel_size
		)

		self$decoder_f <- VplotDecoder(
			vplot_width = self$n_bins_per_block,
			vplot_height = self$n_intervals,
			filters0 = filters0,
			upsample_layers = upsample_layers,
			filters = filters,
			kernel_size = kernel_size
		)

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(latent_dim)),
			scale_identity_multiplier = 1
		)

		function(x, ..., training = TRUE){

			y <- x$vplots %>% 
				self$encoder() 

			y_f <- x$vplots_f %>% 
				self$encoder_f()

			y <- list(y, y_f) %>% tf$concat(axis = 1L)
			y <- y %>% self$dense_1()

      q_m <- y[, 1:self$latent_dim]
			q_v <- tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3)

			posterior <- tfp$distributions$MultivariateNormalDiag(
				loc = q_m,
				scale_diag = q_v
			)

			if (training){
				z <- posterior$sample()
			}else{
				z <- posterior$mean()
			}

			vplots <- z %>% self$decoder(training = training)
			vplots <- vplots %>% tf$keras$activations$softmax(1L)

			vplots_f <- z %>% self$decoder_f(training = training)

			list(
				posterior = posterior, 
				z = z, 
				vplots = vplots,
				vplots_f = vplots_f
			)
		}
	})
}


#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
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
		model = 'FourierVaeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		...
	){

		d <- list()

		vplots <- assays(x)$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
			scale_vplot()

		vplots_f <- vplots %>%
			tf$squeeze(3L) %>%
			tf$cast(dtype = tf$dtypes$complex64) %>%
			tf$signal$fft2d() %>%
			tf$signal$fftshift() %>%
			tf$cast(dtype = tf$dtypes$float32) %>%
			tf$expand_dims(3L) %>%
			tf$image$per_image_standardization()

		list(vplots = vplots, vplots_f = vplots_f)
	}
)


#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @return a VaeModel
#'
setMethod(
	'fit',
	signature(
		model = 'FourierVaeModel',
		x = 'tf_dataset'
	),
	function(
		 model,
		 x,
		 batch_size = 128L,
		 epochs = 100L,
		 learning_rate = 1e-3,
		 compile = TRUE,
		 beta = 5e-5,
		 gamma = 1
	 ){

		if (length(beta) == 1)
			beta <- rep(beta, epochs)
		beta <- tf$cast(beta, tf$float32)

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot
		reconstrution_f_loss <- tf$keras$losses$MeanSquaredError(reduction = 'none')	# loss for the V-plot

		train_step <- function(x, b){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, training = TRUE)
				loss_reconstruction <- (reconstrution_loss(x$vplots, res$vplots)) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_f_reconstruction <- reconstrution_f_loss(x$vplots_f, res$vplots_f) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss_f_reconstruction <- gamma * loss_f_reconstruction
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss_kl <- b * loss_kl
				loss <- loss_reconstruction + loss_kl + loss_f_reconstruction
	 		})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss,
				loss_reconstruction = loss_reconstruction,
				loss_f_reconstruction = loss_f_reconstruction,
				loss_kl = loss_kl
			)
		}


		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		for (epoch in seq_len(epochs)){
			loss <- NULL
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch, beta[epoch])
				loss <- rbind(loss, sapply(res, as.numeric))
			})

			if (epoch == 1 || epoch %% 10 == 0){
				loss <- colMeans(loss)
				sprintf('epoch=%6.d/%6.d | beta=%.3e | %s', epoch, epochs, beta[epoch], paste(sapply(1:length(loss), function(i) sprintf('%s=%13.7f', names(loss)[i], loss[i])), collapse = ' | ')) %>%
					message()
			}
		}
		model
	}
)

#' predict
#'
