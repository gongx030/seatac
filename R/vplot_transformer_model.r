setOldClass('rpytools.call.Transformer')
setOldClass('python.builtin.Transformer')
setClassUnion('rpytools.call.TransformerOrpython.builtin.Transformer', members = c('rpytools.call.Transformer', 'python.builtin.Transformer'))

setClass(
	'vplot_transformer_model',
	slot = c(
		transformer = 'rpytools.call.TransformerOrpython.builtin.Transformer',
		block_size = 'integer'
	)
)



TransformerEncoder <- reticulate::PyClass(
	'TransformerEncoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			num_layers, 
			d_model,
			num_heads, 
			dff, 
			kmers_size, 
			maximum_position_encoding = 2000L,
			rate = 0.1
		){
			super()$`__init__`()

			self$d_model <- d_model
			self$num_layers <- num_layers

			self$embedding <- tf$keras$layers$Embedding(kmers_size, d_model) 

			self$pos_encoding <- positional_encoding(maximum_position_encoding, self$d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)

			self$enc_layers <- lapply(seq_len(num_layers), function(i) TransformerEncoderLayer(d_model, num_heads = num_heads, dff = dff, rate = rate))

			self$dropout <- tf$keras$layers$Dropout(rate)

			NULL

		},
		call = function(self, x, g, training = TRUE, mask = NULL){

			seq_len <- x$shape[2]

			# adding embedding and position encoding.
			x <- self$embedding(x)  # (batch_size, input_seq_len, d_model)
			x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))

			x <- x + self$pos_encoding[, 1:seq_len, ]

			x <- self$dropout(x, training = training)
							    
			for (i in seq_len(self$num_layers)){
				x <- self$enc_layers[[i - 1]](x, training = training, mask = mask)
			}
			x 
		}
	)
)


#' Decoder
#' 
#' A transformer decoder to recover an image
#'
TransformerDecoder <- reticulate::PyClass(
	'TransformerDecoder',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(
			self, 
			num_layers,
			d_model,
			num_heads,
			dff,
			maximum_position_encoding = 2000L,
			rate = 0.1
		){

			super()$`__init__`()

			self$d_model <- d_model
			self$num_layers <- num_layers
			    
			self$pos_encoding <- positional_encoding(maximum_position_encoding, d_model) %>%
				tf$cast(tf$float32) %>%
				tf$expand_dims(0L)
					    
			self$dec_layers <- lapply(seq_len(num_layers), function(i) TransformerDecoderLayer(d_model, num_heads = num_heads, dff = dff, rate = rate) )
			self$dropout <- tf$keras$layers$Dropout(rate)

			NULL
		},
		call = function(self, x, enc_output, training, look_ahead_mask, padding_mask){

			attention_weights <- list()

			seq_len <- x$shape[2]

			x <- x * tf$math$sqrt(tf$cast(self$d_model, tf$float32))
			x <- x + self$pos_encoding[, 1:seq_len, ]
			x <- self$dropout(x, training = training)

			for (i in seq_len(self$num_layers)){
				res <- self$dec_layers[[i - 1]](x, enc_output = enc_output, training = training, look_ahead_mask = look_ahead_mask, padding_mask = padding_mask)
				x <- res$output
				attention_weights[[sprintf('decoder_layer%d_block1', i)]] <- res$attention_weights_block1
				attention_weights[[sprintf('decoder_layer%d_block2', i)]] <- res$attention_weights_block2
			}

			list(
				output = x, 
				attention_weights = attention_weights
			)
		}
	)
)

Transformer <- reticulate::PyClass(
	'Transformer',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self, 
			num_layers,
			d_model,
			num_heads,
			dff,
			kmers_size,
			n_intervals,
			rate = 0.1
		){

			super()$`__init__`()

			self$conv_1 <- tf$keras$layers$Conv2D(
				filters = d_model,
				kernel_size = c(2L, 2L),
				strides = c(2L, 2L),
				activation = 'relu'
			)

			self$conv_2 <- tf$keras$layers$Conv2D(
				filters = d_model,
				kernel_size = c(2L, 2L),
				strides = c(2L, 2L),
				activation = 'relu'
			)

			self$conv_3 <- tf$keras$layers$Conv2D(
				filters = d_model,
				kernel_size = c(2L, 2L),
				strides = c(2L, 2L)
			)

			self$encoder <- TransformerEncoder(
				num_layers = num_layers, 
				d_model = d_model, 
				num_heads = num_heads,  
				dff = dff, 
				kmers_size = kmers_size
			)

			self$decoder <- TransformerDecoder(
				num_layers = num_layers, 
				d_model = d_model, 
				num_heads = num_heads, 
				dff = dff
			)

			self$deconv_1 <- tf$keras$layers$Conv2DTranspose(
				filters = d_model,
				kernel_size = c(2L, 2L),
				strides = c(2L, 2L),
				activation = 'relu'
			)	

			self$deconv_2 <- tf$keras$layers$Conv2DTranspose(
				filters = d_model,
				kernel_size = c(2L, 2L),
				strides = c(2L, 2L),
				activation = 'relu'
			)	

			self$deconv_3 <- tf$keras$layers$Conv2DTranspose(
				filters = 1L,
				kernel_size = c(2L, 2L),
				strides = c(2L, 2L),
				activation = 'relu'
			)	

			self$bn_1 <- tf$keras$layers$BatchNormalization()
			self$bn_2 <- tf$keras$layers$BatchNormalization()
			self$bn_3 <- tf$keras$layers$BatchNormalization()
			self$bn_4 <- tf$keras$layers$BatchNormalization()
			self$bn_5 <- tf$keras$layers$BatchNormalization()

			NULL
		},
		call = function(self, inp, tar, training, enc_padding_mask, dec_padding_mask){

			tar <- tar %>% 
				self$conv_1() %>%
				self$bn_1() %>%
				self$conv_2() %>%
				self$bn_2() %>%
				self$conv_3()

			height <- tar$shape[2]
			width <- tar$shape[3]

			tar <- tar %>% 
				tf$reshape(shape(tar$shape[1], height * width, tar$shape[4]))

#			look_ahead_mask <- create_look_ahead_mask(tar$shape[1])	

			enc_output <- self$encoder(
				inp, 
				training = training, 
				mask = enc_padding_mask
			)  # (batch_size, inp_seq_len, d_model)

			res <- self$decoder(
				tar, 
				enc_output = enc_output, 
				training = training, 
				look_ahead_mask = NULL, 
				padding_mask = dec_padding_mask
			)

			final_output <- res$output %>%
				tf$reshape(shape(res$output$shape[1], height, width, self$encoder$d_model)) %>%
				self$deconv_1() %>%
				self$bn_4() %>%
				self$deconv_2() %>%
				self$bn_5() %>%
				self$deconv_3()

			list(output = final_output, attention_weights = res$attention_weights)

		}
	)
)


#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_transformer_model',
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

		n_bins_per_block <- as.integer(model@block_size / x@bin_size)

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

				xs <- x[h] %>% prepare_blocks_with_kmers(model@block_size + 1L * x@bin_size, min_reads_per_block) 
				n_blocks <- xs$vplot$shape[1]

				for (i in seq_len(steps_per_epoch)){

					b <- sample(0:(n_blocks - 1), batch_size_blocks, replace = TRUE)
					xi <- tf$gather(xs$vplot, b)
					ci <- tf$gather(xs$kmers, b)
					xi_inp <- xi[, , 1:n_bins_per_block, ]
					xi_real <- xi[, , 2:(n_bins_per_block + 1), ]

					with(tf$GradientTape(persistent = TRUE) %as% tape, {

						res <- model@transformer(ci, xi_inp, training = TRUE, enc_padding_mask = NULL, dec_padding_mask = NULL)

						loss <- loss_object(res$output, xi_real) %>%
							tf$reduce_sum(c(1L, 2L)) %>%
							tf$reduce_mean()

					})

					total_loss <- total_loss + loss

					gradients <- tape$gradient(loss, model@transformer$trainable_variables)
					list(gradients, model@transformer$trainable_variables) %>%
						purrr::transpose() %>%
						optimizer$apply_gradients()

				}

				flog.info(sprintf('training %s | epoch=%4.d/%4.d | window batch=%5.d/%5.d | n_blocks=%7.d | total_loss=%13.3f', class(model), epoch, epochs, j, n_batch, n_blocks, total_loss))
			}

			# evaluating the predicted performance
			valid <- sample(which(rowMeans(x$mnase) > 0.5), 500) 
			x_pred <- model %>% predict(x[valid])
			x_pred <- add_nucleosome_signal(x_pred)
			vplot(x_pred, 'predicted_counts', main = sprintf('epoch=%d', epoch))
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
		model = 'vplot_transformer_model',
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

		n_bins_per_block <- as.integer(model@block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		predicted_counts <- matrix(0, length(x), x@n_bins_per_window * x@n_intervals)

		for (i in 1:n_batch){

			if (i == 1 || i %% 100 == 0)
				flog.info(sprintf('predicting %s | batch=%4.d/%4.d', class(model), i, n_batch))

			b <- starts[i]:ends[i]

			inputs <- x[b] %>% prepare_blocks_with_kmers(model@block_size, min_reads = 0L) 
			xi <- inputs$vplots
			ci <- inputs$kmers

			res <- model@transformer(ci, xi, training = TRUE, enc_padding_mask = NULL, dec_padding_mask = NULL)
			xi_pred <- res$output %>%
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

