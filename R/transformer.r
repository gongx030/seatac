#' scaled_dot_product_attention
#'
#' Calculate the attention weights.
#'
#' 	 q, k, v must have matching leading dimensions.
#'   k, v must have matching penultimate dimension, i.e.: seq_len_k = seq_len_v.
#'	 The mask has different shapes depending on its type(padding or look ahead) 
#'	 but it must be broadcastable for addition.
#' 
#' @param q: query shape == (..., seq_len_q, depth)
#' @param k: key shape == (..., seq_len_k, depth)
#' @param v: value shape == (..., seq_len_v, depth_v)
#' @param mask: Float tensor with shape broadcastable to (..., seq_len_q, seq_len_k). Defaults to None.
#'
#' @return output, attention_weights

scaled_dot_product_attention <- function(q, k, v, mask = NULL){

	matmul_qk <- tf$matmul(q, k, transpose_b = TRUE)  # (..., seq_len_q, seq_len_k)
	dk <- tf$cast(tf$shape(k) %>% tail(1), tf$float32)
	scaled_attention_logits <- matmul_qk / tf$math$sqrt(dk)

	if (!is.null(mask)){
		scaled_attention_logits <- scaled_attention_logits + (mask * -1e9)  
	}
	
	attention_weights <- tf$nn$softmax(scaled_attention_logits, axis = -1)  # (..., seq_len_q, seq_len_k)

	output <- tf$matmul(attention_weights, v)  # (..., seq_len_q, depth_v)

	list(output = output, attention_weights = attention_weights)

} # scaled_dot_product_attention

MultiHeadAttention <- reticulate::PyClass(
	'MultiHeadAttention',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(self, d_model, num_heads) {

			super()$`__init__`()

			self$num_heads <- num_heads
			self$d_model <- d_model

			stopifnot(d_model %% self$num_heads == 0)

			self$depth <- as.integer(d_model / self$num_heads)

			self$wv <- tf$keras$layers$Dense(d_model, name = 'v')
			self$wq <- tf$keras$layers$Dense(d_model, name = 'q')
			self$wk <- tf$keras$layers$Dense(d_model, name = 'k')

			self$dense = tf$keras$layers$Dense(d_model, name = 'dense')

			NULL

		},
		split_heads = function(self, x, batch_size){

		  # Split the last dimension into (num_heads, depth).
			#	Transpose the result such that the shape is (batch_size, num_heads, seq_len, depth)
			x <- tf$reshape(x, c(batch_size, -1L, self$num_heads, self$depth))
			tf$transpose(x, perm = c(0L, 2L, 1L, 3L))
		},
		call = function(self, v, k, q, mask = NULL, ...){

			batch_size <- tf$shape(q)[1]
			
			q <- self$wq(q)
			k <- self$wk(k)
			v <- self$wv(v)

			q <- self$split_heads(q, batch_size)  # (batch_size, num_heads, seq_len_q, depth)
			k <- self$split_heads(k, batch_size)  # (batch_size, num_heads, seq_len_q, depth)
			v <- self$split_heads(v, batch_size)  # (batch_size, num_heads, seq_len_q, depth)

			y <- scaled_dot_product_attention(q, k, v, mask)
			scaled_attention <- y$output

			scaled_attention <- tf$transpose(scaled_attention, perm = c(0L, 2L, 1L, 3L))  # (batch_size, seq_len_q, num_heads, depth)
			concat_attention <- tf$reshape(scaled_attention, c(batch_size, -1L, self$d_model))  # (batch_size, seq_len_q, d_model)

			output <- self$dense(concat_attention)  # (batch_size, seq_len_q, d_model)

			list(output = output, attention_weights = y$attention_weights)
		}
	)
)


point_wise_feed_forward_network <- function(d_model, dff){
	model <- tf$keras$Sequential()
	model$add(tf$keras$layers$Dense(dff, activation = 'relu')) # (batch_size, seq_len, dff)
	model$add(tf$keras$layers$Dense(d_model)) # (batch_size, seq_len, d_model)
	model
}


EncoderLayer <- reticulate::PyClass(
	'EncoderLayer',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(self, d_model, num_heads, dff, rate = 0.1) {

			super()$`__init__`()

			self$mha <- MultiHeadAttention(d_model, num_heads)
			self$ffn <- point_wise_feed_forward_network(d_model, dff)

			self$layernorm1 <- tf$keras$layers$LayerNormalization(epsilon=1e-6)
			self$layernorm2 <- tf$keras$layers$LayerNormalization(epsilon=1e-6)
						    
			self$dropout1 <- tf$keras$layers$Dropout(rate)
			self$dropout2 <- tf$keras$layers$Dropout(rate)

			NULL

		},
		call = function(self, x, training, mask){

			res <- self$mha(x, x, x, mask)  # (batch_size, input_seq_len, d_model)
			attn_output <- res$output

			attn_output <- self$dropout1(attn_output, training = training)
			out1 <- self$layernorm1(x + attn_output)  # (batch_size, input_seq_len, d_model)
				    
			ffn_output <- self$ffn(out1)  # (batch_size, input_seq_len, d_model)
			ffn_output = self$dropout2(ffn_output, training = training)
			out2 <- self$layernorm2(out1 + ffn_output)  # (batch_size, input_seq_len, d_model)
			out2
						 
		}
	)
)


Encoder <- reticulate::PyClass(
	'Encoder',
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

			self$enc_layers <- lapply(seq_len(num_layers), function(i) EncoderLayer(d_model, num_heads, dff, rate))

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
Decoder <- reticulate::PyClass(
	'Decoder',
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

