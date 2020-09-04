setOldClass('rpytools.call.CVaeModel')
setOldClass('python.builtin.CVaeModel')
setClassUnion('CVaeModel', members = c('rpytools.call.CVaeModel', 'python.builtin.CVaeModel'))


CVaeEncoder <- reticulate::PyClass(
	'CVaeEncoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			latent_dim,
			n_intervals,
			bin_size,
			kmers_size,
			kmers_kernel_size = 5L,
			rate = 0.1
		){
			super()$`__init__`()

			self$latent_dim <- latent_dim
			self$n_intervals <- as.integer(n_intervals)
			self$bin_size <- as.integer(bin_size)
			self$kmers_size <- as.integer(kmers_size)

			self$dropout_1 <- tf$keras$layers$Dropout(rate)

			self$conv_1d <- tf$keras$layers$Conv1D(
				filters = self$n_intervals, 
				kernel_size = kmers_kernel_size,
				strides = self$bin_size, 
				padding = 'same',
				activation = 'softmax'
			)

			self$vplot_encoder <- VplotEncoder()

			self$dense_1 <- tf$keras$layers$Dense(
				units = 2 * self$latent_dim
			)

			NULL

		},
		call = function(self, x, g, training = TRUE, mask = NULL){

			g <- g %>% 
				tf$one_hot(depth = self$kmers_size) %>%	# one hot encoding of nucleotides or k-mers
				self$conv_1d() %>%
	    	self$dropout_1() %>%
				tf$transpose(shape(0L, 2L, 1L)) %>%
				tf$expand_dims(3L)

			x <- x + g 

			yx <- x %>% 
				self$vplot_encoder() %>%
				self$dense_1()

			yg <- g %>% 
				self$vplot_encoder() %>%
				self$dense_1()

			posterior <- tfp$distributions$MultivariateNormalDiag(
				loc = yx[, 1:self$latent_dim],
				scale_diag = tf$nn$softplus(yx[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3),
				name = 'posterior'
			)

			prior <- tfp$distributions$MultivariateNormalDiag(
				loc = yg[, 1:self$latent_dim],
				scale_diag = tf$nn$softplus(yx[, (self$latent_dim + 1):(2 * self$latent_dim)] + 1e-3),
				name = 'prior'
			)
							    
			list(posterior = posterior, prior = prior)

		}
	)
)


#' Decoder
#' 
#' A transformer decoder to recover an image
#'
CVaeDecoder <- reticulate::PyClass(
	'CVaeDecoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			vplot_width,	# width of the image
			vplot_height, # height of the image
			filters0 = 64,
			rate = 0.1
		){

			super()$`__init__`()

			self$vplot_decoder <- VplotDecoder(
				filters0 = filters0,																								 
				vplot_width = vplot_width, 
				vplot_height = vplot_height
			)

			NULL
		},
		call = function(self, x, training = TRUE){

			y <- x %>%
				self$vplot_decoder() %>%
				tf$keras$activations$softmax(axis = 1L)
			y

		}
	)
)

CVaeModel <- reticulate::PyClass(
	'CVaeModel',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			latent_dim,
			n_intervals,
			bin_size,
			kmers_size,
			block_size,
			filters0 = 64L,
			kmers_kernel_size = 5L,
			rate = 0.1
		){

			super()$`__init__`()

			if (block_size %% bin_size != 0)
				stop('block_size must be a multiple of bin_size')

			self$block_size <- block_size
			self$bin_size <- bin_size
			self$n_bins_per_block <- as.integer(block_size / bin_size)

			self$encoder <- CVaeEncoder(
				latent_dim = latent_dim,
				n_intervals = n_intervals,
				bin_size = bin_size,
				kmers_size = kmers_size,
				kmers_kernel_size = kmers_kernel_size,
				rate = rate
			)

			self$decoder <- CVaeDecoder(
				filters0 = filters0,
				vplot_width = self$n_bins_per_block,
				vplot_height = n_intervals
			)

			NULL
		},
		call = function(self, x, g, training = TRUE){

			res <- self$encoder(x, g)
			z <- res$posterior$sample()
			c <- res$prior$sample()
			x_pred <- (z + c) %>% model$decoder()

			list(
				posterior = res$posterior, 
				prior = res$prior,
				z = z, 
				x_pred = x_pred
			)
		}
	)
)

#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'CVaeModel',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 32L,
		batch_size_predict = 1L,
		batch_size_select = 128L,
		epochs = 100L,
		min_reads_per_block = 10,
		loss_type = 'all',
		plot = FALSE
	){

		optimizer <- tf$keras$optimizers$Adam(learning_rate)

		kl <- tf$keras$losses$KLDivergence(reduction = 'none')

		inputs <- select_blocks(x, block_size = block_size, min_reads = min_reads_per_block, batch_size = batch_size_select, with_kmers = TRUE)

		n <- inputs$vplots$shape[[1]] # total number of samples

		starts <- seq(1, n, by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > n] <- n
		n_batch <- length(starts)

		for (epoch in seq_len(epochs)){

			for (j in seq_len(n_batch)){

				total_loss <- 0
				total_loss_reconstruction <- 0
				total_loss_kl <- 0

				b <- starts[j]:ends[j]

				xi <- tf$gather(inputs$vplot, b)
				xi <- xi %>% scale_vplot()

				ci <- tf$gather(inputs$kmers, b)

				if (loss_type == 'all'){
					w <- tf$ones_like(xi)
				}else if (loss_type == 'non_empty'){
					w <- tf$reduce_sum(xi, 1L, keepdims = TRUE) > 0
					w <- w %>% tf$cast(tf$float32)
				}
				w <- w %>% tf$squeeze(3L)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					res <- model(xi, ci)

					loss_reconstruction <- (w * kl(xi, res$x_pred)) %>%
						tf$reduce_sum(shape(1L, 2L)) %>%
						tf$reduce_mean()

					loss_kl <- (res$posterior$log_prob(res$z) - res$prior$log_prob(res$z)) %>%
						tf$reduce_mean()

					loss <- loss_reconstruction + loss_kl

				})


				total_loss_reconstruction  <- total_loss_reconstruction + loss_reconstruction
				total_loss_kl <- total_loss_kl + loss_kl
				total_loss <- total_loss + loss

				gradients <- tape$gradient(loss, model$trainable_variables)
				list(gradients, model$trainable_variables) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

				flog.info(sprintf('training | epoch=%4.d/%4.d | batch=%5.d/%5.d | total_loss_reconstruction=%13.3f | total_loss_kl=%13.3f | total_loss=%13.3f', epoch, epochs, j, n_batch, total_loss_reconstruction, total_loss_kl, total_loss))

			}


			if (epoch %% 5 == 0){
				# evaluating the predicted performance
				x_sample <- sample(x, 500)
				x_pred <- model %>% predict(x_sample, batch_size = batch_size_predict)
				x_pred <- add_nucleosome_signal(x_pred)
				if (plot){
#					vplot(x_pred, 'counts', main = sprintf('epoch=%d', epoch))
					vplot(x_pred, 'predicted_counts', main = sprintf('epoch=%d', epoch))
				}
				rmse <- sqrt(rowSums((x_pred$full_nucleoatac_scaled - x_pred$nucleosome_signal)^2))
				rmse_nucleoatac <- sqrt(rowSums((x_pred$full_nucleoatac_scaled - x_pred$nucleoatac_scaled)^2))
				flog.info(sprintf('evaluating | epoch=%4.d/%4.d | rmse=%.3f | rmse_nucleoatac=%.3f',  epoch, epochs, mean(rmse), mean(rmse_nucleoatac)))
			}
		}
		model
	}
)


#'
setMethod(
	'predict',
	signature(
		model = 'CVaeModel',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		batch_size = 1L # v-plot per batch
	){

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		n_blocks_per_window <- as.integer(x@n_bins_per_window - model$n_bins_per_block + 1)
		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting | batch=%4.d/%4.d', i, n_batch))

			b <- starts[i]:ends[i]

			inputs <- x[b] %>% prepare_blocks(model$block_size, with_kmers = TRUE)
			xi <- inputs$vplots %>%
				scale_vplot()
			ci <- inputs$kmers

			res <- model(xi, ci)

			xi_pred <- res$x_pred %>%
				tf$reshape(c(length(b), n_blocks_per_window, x@n_intervals, model$n_bins_per_block, 1L)) %>%
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

