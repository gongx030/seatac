setOldClass('rpytools.call.Ac2Model')
setOldClass('python.builtin.Ac2Model')
#setClassUnion('Ac2Model', members = c('rpytools.call.Ac2Model', 'python.builtin.Ac2Model'))
setClassUnion('Ac2Model', members = c('python.builtin.Ac2Model', 'rpytools.call.Ac2Model'))

KmersEncoder <- reticulate::PyClass(
	'KmersEncoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			latent_dim,
			num_layers, 
			d_model, 
			num_heads, 
			dff, 
			kmers_size, 
			bin_size,
			max_pos = 1000L,
			rate = 0.1
		){
			super()$`__init__`()

			self$d_model <- d_model
			self$latent_dim <- latent_dim
			self$num_layers <- num_layers
			self$kmers_size <- kmers_size

			self$embedding <- tf$keras$layers$Embedding(self$kmers_size, d_model)

			self$pool_1 <- tf$keras$layers$MaxPool1D(pool_size = bin_size, strides = bin_size)

			self$pos_encoding <- positional_encoding(max_pos, self$d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)

			self$dropout_1 <- tf$keras$layers$Dropout(rate)

			self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(d_model, num_heads, dff, rate))

			self$dense_1 <- tf$keras$layers$Dense(units = self$latent_dim)

			NULL

		},
		call = function(self, x, training = TRUE, mask = NULL){

			block_size <- x$shape[2]

			# adding embedding and position encoding.
			x <- self$embedding(x)  # (batch_size, block_size, d_model)
			x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))
			x <- x + self$pos_encoding[, 1:block_size, ]
	    x <- self$dropout_1(x, training = training)

			x <- x %>% self$pool_1()

	    for (i in seq_len(self$num_layers)){
	      x <- self$enc_layers[[i - 1]](x, training = training, mask = mask)
			}

			x <- x %>% 
				tf$reduce_mean(2L) %>%
				self$dense_1()
			x

		}
	)
)


Ac2Model<- reticulate::PyClass(
	'Ac2Model',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			latent_dim,
			num_layers,
			d_model,
			num_heads,
			dff,
			kmers_size,
			block_size,
			bin_size,
			n_intervals,
			rate = 0.1
		){

			super()$`__init__`()

			self$block_size <- block_size
			self$bin_size <- bin_size
			self$n_intervals <- n_intervals

			self$vplot_width <- as.integer(block_size / bin_size)
			self$vplot_encoder <- VplotEncoder(latent_dim = latent_dim)

			self$vplot_decoder <- VplotDecoder(
				vplot_width = self$vplot_width,
				vplot_height = self$n_intervals
			)
			
			self$kmers_encoder <- KmersEncoder(
				latent_dim = latent_dim, 
				num_layers = num_layers, 
				d_model = d_model, 
				num_heads = num_heads, 
				dff = dff, 
				kmers_size = length(x@kmers),
				bin_size = bin_size
			)

			NULL
		},
		call = function(self, x, g, training = TRUE, mask = NULL){

			z_x <- self$vplot_encoder(x)
			z_g <- self$kmers_encoder(g)
			z <- z_x + z_g
			y <- self$vplot_decoder(z)
			y
		}
	)
)



#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'Ac2Model',
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
		min_reads_per_block = 10,
		theta = 0
	){

		optimizer <- tf$keras$optimizers$Adam(learning_rate)

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		loss_object <-  tf$keras$losses$MeanSquaredError(reduction = 'none')

		for (epoch in seq_len(epochs)){

			for (j in seq_len(n_batch)){

				total_loss <- 0

				h <- starts[j]:ends[j]

				xs <- x[h] %>% prepare_blocks_with_kmers(model$block_size, min_reads_per_block) 
				n_blocks <- xs$vplot$shape[[1]]

				for (i in seq_len(steps_per_epoch)){

					b <- sample(0:(n_blocks - 1), batch_size_blocks, replace = TRUE)
					xi <- tf$gather(xs$vplot, b)
					ci <- tf$gather(xs$kmers, b)

					mask <- tf$random$uniform(xi$shape) > theta
					mask <- mask %>% tf$cast(tf$float32)
					xi_input <- xi * mask

					with(tf$GradientTape(persistent = TRUE) %as% tape, {

						xi_pred <- model(xi_input, ci)

						loss <- loss_object(xi, xi_pred) %>%
							tf$reduce_sum(c(1L, 2L)) %>%
							tf$reduce_mean()

					})

					total_loss <- total_loss + loss

					gradients <- tape$gradient(loss, model$trainable_variables)
					list(gradients, model$trainable_variables) %>%
						purrr::transpose() %>%
						optimizer$apply_gradients()

				}

				flog.info(sprintf('training | theta=%.3f | epoch=%4.d/%4.d | window batch=%5.d/%5.d | n_blocks=%7.d | total_loss=%13.3f', theta, epoch, epochs, j, n_batch, n_blocks, total_loss))
			}

			# evaluating the predicted performance
			valid <- sample(which(rowMeans(x$mnase) > 0.5), 500)
			x_pred <- model %>% predict(x[valid])
			x_pred <- add_nucleosome_signal(x_pred)
#			vplot(x_pred, 'predicted_counts', main = sprintf('epoch=%d', epoch))
			rmse <- sqrt(rowSums((x_pred$mnase_scaled - x_pred$nucleosome_signal)^2))
			flog.info(sprintf('evaluating | epoch=%4.d/%4.d | rmse=%.3f',  epoch, epochs, mean(rmse)))

		}
		model
	}
)


#'
setMethod(
	'predict',
	signature(
		model = 'Ac2Model',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		batch_size = 4L # v-plot per batch
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		block_size <- model$block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting | batch=%4.d/%4.d', i, n_batch))

			b <- starts[i]:ends[i]

			inputs <- x[b] %>% prepare_blocks_with_kmers(model$block_size)
			xi <- inputs$vplots
			ci <- inputs$kmers

			xi_pred <- model(xi, ci) %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, n_bins_per_block, 1L)) %>%
				reconstruct_vplot_from_blocks() %>%
				tf$squeeze(-1L) %>%
				as.array()

			xi_pred <- aperm(xi_pred, c(1, 3, 2))
			dim(xi_pred) <- c(length(b), x@n_bins_per_window * x@n_intervals)

			predicted_counts[b, ] <- xi_pred

		}

		mcols(x)$predicted_counts <- predicted_counts

		x
	}
) # predict

