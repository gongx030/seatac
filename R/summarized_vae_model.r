
#' SummarizedVaeModel
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
SummarizedVaeModel <- function(
	 latent_dim = 10L,
	 n_intervals = 48L,
	 bin_size = 5L,
	 block_size = 640L,
	 filters0 = 128L,
	 encoder_filters = c(32L, 32L, 32L),
	 encoder_kernel_size = c(3L, 3L, 3L),
	 encoder_window_strides = c(2L, 2L, 2L),
	 encoder_interval_strides = c(2L, 2L, 2L),
	 decoder_filters = c(32L, 32L, 1L),
	 decoder_kernel_size = c(3L, 3L, 3L),
	 decoder_window_strides = c(2L, 2L, 2L),
	 decoder_interval_strides = c(2L, 2L, 1L),
	 rate = 0.1,
	 name = NULL
){

	keras_model_custom(name = name, function(self){

		if (block_size %% bin_size != 0)
			stop('block_size must be a multiple of bin_size')

		self$block_size <- block_size
		self$bin_size <- bin_size
		self$n_bins_per_block <- as.integer(block_size / bin_size)

		self$encoder <- VaeEncoder(
			latent_dim = latent_dim,
			filters = encoder_filters,
			kernel_size = encoder_kernel_size,
			window_strides = encoder_window_strides,
			interval_strides = encoder_interval_strides,
			distribution = 'MultivariateNormalDiag'
		)

		self$library_encoder <- VaeEncoder(
			latent_dim = 1L,
			filters = 32L,
			kernel_size = 3L,
			window_strides = 2L,
			interval_strides = 2L,
			distribution = 'MultivariateNormalDiag'
		)

		self$decoder <- VplotDecoder(
			vplot_width = self$n_bins_per_block,
			vplot_height = n_intervals,
			filters0 = filters0,
			filters = decoder_filters,
			kernel_size = decoder_kernel_size,
			window_strides = decoder_window_strides,
			interval_strides = decoder_interval_strides
		)

		self$prior <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(latent_dim)),
			scale_identity_multiplier = 1
		)

		self$prior_library <- tfp$distributions$MultivariateNormalDiag(
			loc  = tf$zeros(shape(1L)),
			scale_identity_multiplier = 1
		)

		function(x, training = TRUE){

			posterior <- self$encoder(x)
			z <- posterior$sample()
			x_pred <- z %>% self$decoder()

			posterior_library <- self$library_encoder(x)
			library_size <- posterior_library$sample()

			lib <- library_size %>%
				tf$expand_dims(2L) %>%
				tf$expand_dims(3L)
			
			list(
				posterior = posterior, 
				posterior_library = posterior_library,
				library_size = library_size,
				z = z, 
				vplots = tfp$distributions$Poisson(log_rate = x_pred + lib)
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

		decoder_filters <- model@model$decoder$filters # a ListWrapper, zero-based
		stopifnot(channels == decoder_filters[[length(decoder_filters) - 1]])	

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
		x = 'tensorflow.python.framework.sparse_tensor.SparseTensor'
	),
	function(
		 model,
		 x,
		 batch_size = 32L,
		 epochs = 100L,
		 learning_rate = 1e-4
	 ){

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		train_step <- function(x){

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x)
				loss_reconstruction <- -res$vplots$log_prob(x) %>%
					tf$reduce_sum(shape(1L, 2L, 3L)) %>%
					tf$reduce_mean()
				loss_kl <- (res$posterior$log_prob(res$z) - model@model$prior$log_prob(res$z)) %>%
					tf$reduce_mean()
				loss_kl_library <- (res$posterior_library$log_prob(res$library_size) - model@model$prior_library$log_prob(res$library_size)) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction + loss_kl + loss_kl_library
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

#		train_step <- tf_function(train_step) # convert to graph mode

		n <- x$shape[[1]]
		n_intervals <- x$shape[[2]]
		n_bins_per_window <- x$shape[[3]]
		channels <- x$shape[[4]]

		for (epoch in seq_len(epochs)){

			batches <- cut_data(n, batch_size)

			# training 
			loss <- NULL
			loss_reconstruction <- NULL
			loss_kl <- NULL

			for (i in 1:length(batches)){

				b <- batches[[i]]

			 	xb <- tf$sparse$slice(x, c(b[1] - 1L, 0L, 0L, 0L),c(length(b), n_intervals, n_bins_per_window, channels))
				xb <- xb %>% tf$sparse$to_dense() 	# to dense

				res <- train_step(xb)

				loss <- c(loss, as.numeric(res$loss))
				loss_reconstruction <- c(loss_reconstruction, as.numeric(res$loss_reconstruction))
				loss_kl <- c(loss_kl, as.numeric(res$loss_kl))
			}

			message(sprintf('epoch=%6.d/%6.d | recon_loss=%15.7f | kl_loss=%15.7f | loss=%15.7f', epoch, epochs, mean(loss_reconstruction), mean(loss_kl), mean(loss)))

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
		x = 'tensorflow.python.framework.sparse_tensor.SparseTensor'
	),
	function(
		model,
		x,
		batch_size = 8L # v-plot per batch
	){

		n <- x$shape[[1]]
		n_intervals <- x$shape[[2]]
		n_bins_per_window <- x$shape[[3]]
		channels <- x$shape[[4]]
		batches <- cut_data(n, batch_size)
		latent <- list()
		library_size <- list()
		for (i in 1:length(batches)){
			b <- batches[[i]]
		 	xb <- tf$sparse$slice(x, c(b[1] - 1L, 0L, 0L, 0L),c(length(b), n_intervals, n_bins_per_window, channels))
			xb <- xb %>% tf$sparse$to_dense() 	# to dense
			posterior <- model@model$encoder(xb)
			library_size[[i]] <- model@model$library_encoder(xb)
			z <- posterior$mean()
			latent[[i]] <- z

		}
		latent <- tf$concat(latent, axis = 0L)
		library_size <- tf$concat(library_size, axis = 0L)
		list(
			latent = latent %>% as.matrix(),
			library_size = library_size %>% as.matrix()
		)
	}
)

