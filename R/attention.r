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
#'
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


#' MultiHeadAttention
#'
MultiHeadAttention <- function(
	d_model, 
	num_heads,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$num_heads <- num_heads
		self$d_model <- d_model

		stopifnot(self$d_model %% self$num_heads == 0)

		self$depth <- as.integer(self$d_model / self$num_heads)
		self$wv <- tf$keras$layers$Dense(self$d_model)
		self$wq <- tf$keras$layers$Dense(self$d_model)
		self$wk <- tf$keras$layers$Dense(self$d_model)
		self$dense <- tf$keras$layers$Dense(self$d_model)

		self$split_heads <- function(x, batch_size){
			# Split the last dimension into (num_heads, depth).
			#	Transpose the result such that the shape is (batch_size, num_heads, seq_len, depth)
			x <- tf$reshape(x, c(batch_size, -1L, self$num_heads, self$depth))
			tf$transpose(x, perm = c(0L, 2L, 1L, 3L))
		}

		function(inputs, mask = NULL, ...){
			batch_size <- tf$shape(inputs$q)[1]
			q <- self$wq(inputs$q)
			k <- self$wk(inputs$k)
			v <- self$wv(inputs$v)
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
	})
}

#' point_wise_feed_forward_network
#'
point_wise_feed_forward_network <- function(d_model, dff){
	model <- tf$keras$Sequential()
	model$add(tf$keras$layers$Dense(dff, activation = 'relu')) # (batch_size, seq_len, dff)
	model$add(tf$keras$layers$Dense(d_model)) # (batch_size, seq_len, d_model)
	model
}
