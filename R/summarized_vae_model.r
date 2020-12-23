#' SummarizedVaeModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
SummarizedVaeModel <- function(
	 latent_dim_fragment_size = 2L,
	 latent_dim_channel = 10L,
	 n_intervals = 48L,
	 bin_size = 5L,
	 block_size = 640L,
	 filters0 = 128L,
	 channels = 1L,
	 encoder_filters = c(32L, 32L),
	 encoder_kernel_size = c(2L, 2L),
	 encoder_window_strides = c(2L, 2L),
	 encoder_interval_strides = c(2L, 2L),
	 decoder_kernel_size = 2L,
	 decoder_strides = 2L,
	 rate = 0.1,
	 name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$block_size <- block_size
		self$bin_size <- bin_size
		self$n_bins_per_block <- as.integer(block_size / bin_size)
		self$channels <- channels
		self$n_intervals <- n_intervals
		self$latent_dim_channel <- latent_dim_channel
		self$latent_dim_fragment_size <- latent_dim_fragment_size
		self$latent_dim <- latent_dim_channel + latent_dim_fragment_size

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(self$latent_dim)),
			scale_identity_multiplier = 1
		)

		self$encoder <- VaeEncoder(
			latent_dim = self$latent_dim,
			filters = encoder_filters,
			kernel_size = encoder_kernel_size,
			window_strides = encoder_window_strides,
			interval_strides = encoder_interval_strides,
			rate = rate,
			distribution = 'MultivariateNormalDiag'
		)

		self$channel_decoder <- VplotDecoder(
			vplot_width = self$n_bins_per_block,
			vplot_height = 1L,
			filters0 = filters0,
			filters = self$channels,
			kernel_size = decoder_kernel_size,
			window_strides = decoder_strides,
			interval_strides = 1L
		)

		self$fragment_size_decoder <- VplotDecoder(
			vplot_width = 1L,
			vplot_height = self$n_intervals,
			filters0 = filters0,
			filters = 1L,
			kernel_size = decoder_kernel_size,
			window_strides = 1L,
			interval_strides = decoder_strides
		)

		function(x, ..., training = TRUE){

			posterior <- self$encoder(x)

			if (training){
				z <- posterior$sample()
			}else{
				z <- posterior$mean()
			}

			z_channel <- z[, 1:self$latent_dim_channel]
			z_fragment_size <- z[, (self$latent_dim_channel + 1):(self$latent_dim_channel + latent_dim_fragment_size)]

			x_pred <- z_channel %>% self$channel_decoder()
			fragment_size <- z_fragment_size %>% self$fragment_size_decoder()

			list(
				posterior = posterior, 
				z = z,
				z_channel = z_channel,
				z_fragment_size = z_fragment_size,
				x_pred = x_pred,
				fragment_size = fragment_size,
				vplots = fragment_size + x_pred
			)
		}
	})
}



#' prepare_data
#'
#' Prepare tfdataset for training and testing a VaeModel model
#'
#' @export
#'
setMethod(
	'prepare_data',
	signature(
		model = 'SummarizedVaeModel',
		x = 'SummarizedVplotsList'
	),
	function(
		model,
		x
	){

		channels <- length(x[[1]])
		n_bins_per_window <- x[[1]]@n_bins_per_window
		n_intervals <- x[[1]]@n_intervals
		n <- length(x)

		x <- do.call('rbind', lapply(x, function(xx) assays(xx)$counts))

		b <- rep(1:n, each = channels)	# row-wise batch index
		m <- rep(1:channels, n)	# row-wise motif index
		w <- rep(1:n_bins_per_window, n_intervals)
		h <- rep(1:n_intervals, each = n_bins_per_window)

		x <- summary(x) %>% as.matrix()

		x <- tf$SparseTensor(
			indices = cbind(b[x[, 1]], m[x[, 1]], w[x[, 2]], h[x[, 2]]) - 1L,
      values = x[, 3],
      dense_shape = c(n, channels, n_bins_per_window, n_intervals)
    ) 
		x <- x %>% tf$sparse$transpose(shape(0L, 3L, 2L, 1L))

		x <- x %>% tf$sparse$to_dense()
		total <- x %>% tf$reduce_sum()
		bs <- x %>% tf$reduce_sum(shape(1L, 2L, 3L), keepdims = TRUE)
		vs <- x %>% tf$reduce_sum(shape(0L, 2L, 3L), keepdims = TRUE)
		hs <- x %>% tf$reduce_sum(shape(0L, 1L, 3L), keepdims = TRUE)
		cs <- x %>% tf$reduce_sum(shape(0L, 1L, 2L), keepdims = TRUE)

		xe <- bs %>%
			tf$math$divide(total) %>%
			tf$multiply(vs) %>%
			tf$math$divide(total) %>%
			tf$multiply(hs) %>%
			tf$math$divide(total) %>%
			tf$multiply(cs)
		x <- (x - xe) / tf$math$maximum(xe, 1)
		x
	}
)


#' fit
#'
#' @export
#'
setMethod(
	'fit',
	signature(
		model = 'SummarizedVaeModel',
		x = 'tf_dataset'
	),
	function(
		 model,
		 x,
		 batch_size = 32L,
		 epochs = 100L,
		 learning_rate = 1e-3,
		 compile = FALSE
	 ){

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		mse <- tf$keras$losses$MeanSquaredError(reduction = 'none')

		x <- x %>%
	    dataset_batch(batch_size)

		train_step <- function(x){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x)
				loss_reconstruction <- mse(tf$expand_dims(x, 4L), tf$expand_dims(res$vplots, 4L)) %>%
					tf$reduce_sum(shape(1L, 2L, 3L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction + loss_kl 
	 		})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss_reconstruction = loss_reconstruction,
				loss_kl = loss_kl,
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

    for (epoch in seq_len(epochs)){

			loss <- NULL
			loss_reconstruction <- NULL
			loss_kl <- NULL
			
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()

			until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch)
				loss <- c(loss, as.numeric(res$loss))
				loss_reconstruction <- c(loss_reconstruction, as.numeric(res$loss_reconstruction))
				loss_kl <- c(loss_kl, as.numeric(res$loss_kl))
			})

			sprintf('epoch=%6.d/%6.d | recon_loss=%15.7f | kl_loss=%15.7f | loss=%15.7f', epoch, epochs, mean(loss_reconstruction), mean(loss_kl), mean(loss)) %>% 
				message()
		}
		model
	}
)

#'
#' @export 
#'
setMethod(
	'predict',
	signature(
		model = 'SummarizedVaeModel',
		x = 'tf_dataset'
	),
	function(
		model,
		x,
		batch_size = 8L # v-plot per batch
	){

		res <- list(
			z = list(),
			z_channel = list(),
			z_fragment_size = list(),
			x_pred = list(),
			vplots = list(),
			fragment_size = list(),
			x_scaled = list()
		)

		iter <- x %>%
	    dataset_batch(batch_size) %>%
			make_iterator_one_shot()

		until_out_of_range({
			batch <- iterator_get_next(iter)
			y <- model@model(batch, training = FALSE)
			y$x_scaled <- batch
			for (j in names(res)){
				res[[j]] <- c(res[[j]], y[[j]])
			}
		})

		for (j in names(res)){
			res[[j]] <- tf$concat(res[[j]], axis = 0L)
		}
		res
	}
)

