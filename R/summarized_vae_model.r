#' SummarizedVplotEncoder
#'
SummarizedVaeEncoder <- function(
	latent_dim = 1L,
	filters = c(32L, 32L, 32L),
	kernel_size = c(3L, 3L, 3L),
	window_strides = c(2L, 2L, 2L),
	interval_strides = c(2L, 2L, 2L),
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self) {

		self$latent_dim <- latent_dim

		self$vplot_encoder <- VplotEncoder(
			filters = filters,
			kernel_size = kernel_size,
			window_strides = window_strides,
			interval_strides = interval_strides
		)

		self$dense_1 <- tf$keras$layers$Dense(1L, activation = 'relu')
    self$dropout_1 <- tf$keras$layers$Dropout(0.8)
		self$dense_2 <- tf$keras$layers$Dense(units = 2 * self$latent_dim)

		function(x, mask = NULL){

			y <- x %>%
				self$vplot_encoder() %>%
				tf$transpose(shape(0L, 3L, 1L, 2L)) %>%
				tf$reshape(shape(x$shape[[1]], self$vplot_encoder$filters[[self$vplot_encoder$n_layers - 1]], -1L)) %>%
				self$dense_1() %>%
				self$dropout_1()  %>%
				tf$squeeze(2L) %>%
				self$dense_2()

			q_m <- y[, 1:self$latent_dim]
			q_v <- tf$nn$softplus(y[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3)
			tfp$distributions$MultivariateNormalDiag(
				loc = q_m,
				scale_diag = q_v
			)
		}
	})
}


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
		self$channels <- decoder_filters[length(decoder_filters)]

		self$encoder <- SummarizedVaeEncoder(
			latent_dim = latent_dim,
			filters = encoder_filters,
			kernel_size = encoder_kernel_size,
			window_strides = encoder_window_strides,
			interval_strides = encoder_interval_strides
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

		self$fragment_size_dense_1 <- tf$keras$layers$Dense(4L, activation = 'relu')
		self$fragment_size_dropout_1 <- tf$keras$layers$Dropout(rate)
		self$fragment_size_dense_2 <- tf$keras$layers$Dense(n_intervals)

		function(x, fragment_size, training = TRUE){

			posterior <- self$encoder(x)
			z <- posterior$sample()
			x_pred <- z %>% self$decoder()

			fragment_size <- fragment_size %>% 
				self$fragment_size_dense_1() %>%
				self$fragment_size_dropout_1() %>%
				self$fragment_size_dense_2() %>%
				tf$expand_dims(2L) %>%
				tf$expand_dims(3L)

			list(
				posterior = posterior, 
				z = z, 
				vplots = x_pred + fragment_size
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
		 learning_rate = 1e-4,
		 compile = FALSE
	 ){

		optimizer <- tf$keras$optimizers$Adam(learning_rate, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)

		reconstrution_loss <- tf$keras$losses$MeanSquaredError(reduction = 'none')	# loss for the V-plot
		n <- x$shape[[1]]

		train_step <- function(x, xn){

			fragment_size <- x %>% 
				tf$reduce_sum(shape(2L, 3L)) %>% 
				scale01()

			xe <- x %>% 
				tf$reduce_sum(shape(1L, 2L, 3L), keepdims = TRUE)	%>%# batch wise read count sum
				tf$multiply(xn) 

			x <- x - xe
			x <- x %>% tf$image$per_image_standardization()

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, fragment_size)
				loss_reconstruction <- reconstrution_loss(x, res$vplots) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
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
				loss = loss,
				loss_reconstruction = loss_reconstruction,
				loss_kl = loss_kl
			)
		}

		if (compile)
			train_step <- tf_function(train_step) # convert to graph mode

		n <- x$shape[[1]]
		n_intervals <- x$shape[[2]]
		n_bins_per_window <- x$shape[[3]]
		channels <- x$shape[[4]]

		total <- x %>% tf$sparse$reduce_sum()
		bs <- x %>% tf$sparse$reduce_sum(shape(1L, 2L, 3L), keepdims = TRUE)	# batch wise read count sum
		vs <- x %>% tf$sparse$reduce_sum(shape(0L, 2L, 3L), keepdims = TRUE)	# vertical (fragment size wise) read count sum
		hs <- x %>% tf$sparse$reduce_sum(shape(0L, 1L, 3L), keepdims = TRUE)	# horizontal (position wise) read count sum
		cs <- x %>% tf$sparse$reduce_sum(shape(0L, 1L, 2L), keepdims = TRUE)	# channel read count sum
		xn <- vs %>%
			tf$math$divide(total) %>%
			tf$multiply(hs) %>% 
			tf$math$divide(total) %>% 
			tf$multiply(cs) %>%
			tf$math$divide(total) 

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
				res <- train_step(xb, xn)

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
		x_pred <- list()
		x_scaled <- list()

		total <- x %>% tf$sparse$reduce_sum()
		bs <- x %>% tf$sparse$reduce_sum(shape(1L, 2L, 3L), keepdims = TRUE)	# batch wise read count sum
		vs <- x %>% tf$sparse$reduce_sum(shape(0L, 2L, 3L), keepdims = TRUE)	# vertical (fragment size wise) read count sum
		hs <- x %>% tf$sparse$reduce_sum(shape(0L, 1L, 3L), keepdims = TRUE)	# horizontal (position wise) read count sum
		cs <- x %>% tf$sparse$reduce_sum(shape(0L, 1L, 2L), keepdims = TRUE)	# channel read count sum

		xn <- vs %>%
			tf$math$divide(total) %>%
			tf$multiply(hs) %>% 
			tf$math$divide(total) %>% 
			tf$multiply(cs) %>%
			tf$math$divide(total) 

#		xn <- vs %>%
#			tf$math$divide(total) %>%
#			tf$multiply(cs) %>%
#			tf$math$divide(total) 

		for (i in 1:length(batches)){

			b <- batches[[i]]
		 	xb <- tf$sparse$slice(x, c(b[1] - 1L, 0L, 0L, 0L),c(length(b), n_intervals, n_bins_per_window, channels))
			xb <- xb %>% tf$sparse$to_dense() 	# to dense

			fragment_size <- xb %>% 
				tf$reduce_sum(shape(2L, 3L), keepdims) %>% 
				scale01()

			xe <- xb %>% 
				tf$reduce_sum(shape(1L, 2L, 3L), keepdims = TRUE)	%>%# batch wise read count sum
				tf$multiply(xn) 

			xb <- xb - xe
			xb <- xb %>% tf$image$per_image_standardization()

			posterior <- model@model$encoder(xb, fragment_size)
			z <- posterior$mean()
			latent[[i]] <- z
			x_pred[[i]] <- model@model$decoder(z)
			x_scaled[[i]] <- xb

		}
		latent <- tf$concat(latent, axis = 0L)
		x_pred <- tf$concat(x_pred, axis = 0L)
		x_scaled <- tf$concat(x_scaled, axis = 0L)

		list(
			latent = latent, 
			x_pred = x_pred,
			x_scaled = x_scaled
		)
	}
)

