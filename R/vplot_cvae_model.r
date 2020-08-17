setOldClass('rpytools.call.Encoder')
setOldClass('python.builtin.Encoder')
setClassUnion('rpytools.call.EncoderOrpython.builtin.Encoder', members = c('rpytools.call.Encoder', 'python.builtin.Encoder'))

setOldClass('rpytools.call.Decoder')
setOldClass('python.builtin.Decoder')
setClassUnion('rpytools.call.DecoderOrpython.builtin.Decoder', members = c('rpytools.call.Decoder', 'python.builtin.Decoder'))

setClass(
	'vplot_cvae_model',
	slot = c(
		encoder = 'rpytools.call.EncoderOrpython.builtin.Encoder',
		decoder = 'rpytools.call.DecoderOrpython.builtin.Decoder'
	)
)



CvaeEncoder <- reticulate::PyClass(
	'CvaeEncoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			latent_dim,
			num_layers, 
			d_model, 
			num_heads, 
			dff, 
			bin_size,
			kmers_size, 
			block_size, 
			rate = 0.1
		){
			super()$`__init__`()

			self$d_model <- d_model
			self$latent_dim <- latent_dim
			self$num_layers <- num_layers
			self$kmers_size <- kmers_size
			self$block_size <- block_size
			self$bin_size <- bin_size

			self$embedding <- tf$keras$layers$Embedding(self$kmers_size, d_model)

			self$pos_encoding <- positional_encoding(self$block_size, self$d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)

			self$dropout_1 <- tf$keras$layers$Dropout(rate)

			self$average_pooling_1d_1<- tf$keras$layers$AveragePooling1D(pool_size = bin_size)

			self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(d_model, num_heads, dff, rate))

			self$conv1d <- tf$keras$layers$Conv1D(
				filters = self$d_model, 
				kernel_size = 3L,
				strides = 1L,
				padding = 'same'
			)

			self$dense_1 <- tf$keras$layers$Dense(
				units = 2 * self$latent_dim
			)

			NULL

		},
		call = function(self, x, g, training = TRUE, mask = NULL){

			# adding embedding and position encoding.
			g <- self$embedding(g)  # (batch_size, block_size, d_model)
			g <- g +  tf$math$sqrt(tf$cast(self$d_model, tf$float32))
			g <- g + self$pos_encoding
	    g <- self$dropout_1(g, training = training)
			g <- g %>% self$average_pooling_1d_1()

			x <- x %>% 
				tf$transpose(c(0L, 2L, 1L, 3L)) %>%
				tf$squeeze(3L) %>%
				self$conv1d()

			x <- x + g 
							    
	    for (i in seq_len(self$num_layers)){
	      g <- self$enc_layers[[i - 1]](g, training = training, mask = mask)
	      x <- self$enc_layers[[i - 1]](x, training = training, mask = mask)
			}

			x <- x %>% tf$reduce_mean(1L)
			g <- g %>% tf$reduce_mean(1L)

			yx <- x %>%
				self$dense_1() 

			yc <- g %>%
				self$dense_1() 

			posterior <- tfp$distributions$MultivariateNormalDiag(
				loc = yx[, 1:self$latent_dim],
				scale_diag = tf$nn$softplus(yx[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3),
				name = 'posterior'
			)

			prior <- tfp$distributions$MultivariateNormalDiag(
				loc = yc[, 1:self$latent_dim],
				scale_diag = tf$nn$softplus(yc[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3),
				name = 'prior'
			)

			list(posterior = posterior, prior = prior, context = g)

		}
	)
)


#' Decoder
#' 
#' A transformer decoder to recover an image
#'
CvaeDecoder <- reticulate::PyClass(
	'CvaeDecoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			vplot_width,	# width of the image
			vplot_height, # height of the image
			rate = 0.1
		){

			super()$`__init__`()

			self$vplot_width <- vplot_width
			self$vplot_height <- vplot_height

			filters0 <- 64L
			filters <- c(32L, 32L, 1L)
			kernel_size <- c(3L, 3L, 3L)
			window_strides <- c(2L, 2L, 2L)
			interval_strides <- c(2L, 2L, 2L)
			window_dim0 <- as.integer(vplot_width / prod(window_strides))
			interval_dim0 <- as.integer(vplot_height / prod(interval_strides))
			output_dim0 <- as.integer(window_dim0 * interval_dim0 * filters0)

			self$dense_1 <- tf$keras$layers$Dense(
				units = output_dim0,
				activation = 'relu'
			)

			self$deconv_1 <- tf$keras$layers$Conv2DTranspose(
				filters = filters[1],
				kernel_size = kernel_size[1],
				strides = shape(interval_strides[1], window_strides[1]),
				padding = 'same',
				activation = 'relu'
			)

			self$deconv_2 <- tf$keras$layers$Conv2DTranspose(
				filters = filters[2],
				kernel_size = kernel_size[2],
				strides = shape(interval_strides[2], window_strides[2]),
				padding = 'same',
				activation = 'relu'
			)

			self$deconv_3 <- tf$keras$layers$Conv2DTranspose(
				filters = filters[3],
				kernel_size = kernel_size[3],
				strides = shape(interval_strides[3], window_strides[3]),
				padding = 'same'
			)

			self$bn_1 <- tf$keras$layers$BatchNormalization()
			self$bn_2 <- tf$keras$layers$BatchNormalization()

			self$reshape_1 <- tf$keras$layers$Reshape(target_shape = c(interval_dim0, window_dim0, filters0))

			NULL
		},
		call = function(self, z, g, training = TRUE){


			y <- tf$concat(list(z, g), axis = 1L) %>%
				self$dense_1() %>%
				self$reshape_1() %>%
				self$deconv_1() %>%
				self$bn_1() %>%
				self$deconv_2() %>%
				self$bn_2() %>%
				self$deconv_3() 

			tfp$distributions$Independent(
				tfp$distributions$Poisson(log_rate = y),
				reinterpreted_batch_ndims = 3L
			)
		}
	)
)

#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_cvae_model',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 32L,
		batch_size_blocks = 256L,
		steps_per_epoch = 16L,
		epochs = 100L,
		min_reads_per_block = 10
	){

		optimizer <- tf$keras$optimizers$Adam(learning_rate)

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		for (epoch in seq_len(epochs)){

			for (j in seq_len(n_batch)){

				total_loss <- 0
				total_loss_reconstruction <- 0
				total_loss_kl <- 0

				h <- starts[j]:ends[j]

				xs <- x[h] %>% prepare_blocks(model, min_reads_per_block) 
				n_blocks <- xs$vplot$shape[0]

				for (i in seq_len(steps_per_epoch)){

					b <- sample(0:(n_blocks - 1), batch_size_blocks, replace = TRUE)
					xi <- tf$gather(xs$vplot, b)
					ci <- tf$gather(xs$kmers, b)

					with(tf$GradientTape(persistent = TRUE) %as% tape, {

						enc <- xi %>% model@encoder(ci)
						z <- enc$posterior$sample()

						xi_pred <- z %>%  
							model@decoder(enc$context)

						loss_kl <- (enc$posterior$log_prob(z) - enc$prior$log_prob(z)) %>%
							tf$reduce_mean()

						loss_reconstruction <- -xi_pred$log_prob(xi) %>%
							tf$reduce_mean()

						loss <- loss_reconstruction + loss_kl

					})

					total_loss_reconstruction  <- total_loss_reconstruction + loss_reconstruction
					total_loss_kl <- total_loss_kl + loss_kl
					total_loss <- total_loss + loss

					encoder_gradients <- tape$gradient(loss, model@encoder$trainable_variables)
					list(encoder_gradients, model@encoder$trainable_variables) %>%
						purrr::transpose() %>%
						optimizer$apply_gradients()

					decoder_gradients <- tape$gradient(loss, model@decoder$trainable_variables)
					list(decoder_gradients, model@decoder$trainable_variables) %>%
						purrr::transpose() %>%
						optimizer$apply_gradients()

				}

				flog.info(sprintf('training %s | epoch=%4.d/%4.d | window batch=%5.d/%5.d | n_blocks=%7.d | total_loss_reconstruction=%13.3f | total_loss_kl=%13.3f | total_loss=%13.3f', class(model), epoch, epochs, j, n_batch, n_blocks, total_loss_reconstruction, total_loss_kl, total_loss))
			}

			# evaluating the predicted performance
			valid <- sample(which(rowMeans(x$mnase) > 0.5 & rowSums(x$counts) < 25), 50)
			x_pred <- model %>% predict(x[valid])
			x_pred <- add_nucleosome_signal(x_pred)
			rmse <- sqrt(rowSums((x_pred$mnase_scaled - x_pred$nucleosome_signal)^2))
			flog.info(sprintf('evaluating %s | epoch=%4.d/%4.d | rmse=%.3f',  class(model), epoch, epochs, mean(rmse)))

		}
		model
	}
)


#'
setMethod(
	'predict',
	signature(
		model = 'vplot_cvae_model',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		batch_size = 8L# v-plot per batch
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		block_size <- model@encoder$block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		latent <- array(NA, c(length(x),  n_blocks_per_window, model@encoder$latent_dim))
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting %s | batch=%4.d/%4.d', class(model), i, n_batch))

			b <- starts[i]:ends[i]

			inputs <- x[b] %>% prepare_blocks(model, min_reads = 0L)
			xi <- inputs$vplots
			ci <- inputs$kmers

			enc <- xi %>% model@encoder(ci)

			z <- enc$posterior$mean() 

			xi_pred <- z %>% model@decoder(enc$context)

			xi_pred <- xi_pred$mean() %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			xi_pred <- aperm(xi_pred, c(1, 3, 2))
			dim(xi_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)

			z <- z %>%
				tf$reshape(c(length(b), n_blocks_per_window, model@encoder$latent_dim))

			latent[b, , ] <- z %>%
				as.array()

			predicted_counts[b, ] <- xi_pred

		}

		mcols(x)$latent <- latent
		mcols(x)$predicted_counts <- predicted_counts

		x
	}
) # predict

