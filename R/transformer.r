#' TransformerEncoderLayer
#'
TransformerEncoderLayer <- function(
	d_model, 
	num_heads, 
	dff, 
	rate = 0.1,
	name = NULL
){
	keras_model_custom(name = name, function(self){
		self$mha <- MultiHeadAttention(d_model, num_heads)
		self$ffn <- point_wise_feed_forward_network(d_model, dff)
		self$layernorm1 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
		self$layernorm2 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
		self$dropout1 <- tf$keras$layers$Dropout(rate)
		self$dropout2 <- tf$keras$layers$Dropout(rate)

		function(x, training = TRUE, mask = NULL){
			res <- self$mha(list(q = x, k = x, v = x), mask)  # (batch_size, input_seq_len, d_model)
			attn_output <- res$output
			attn_output <- self$dropout1(attn_output, training = training)
			out1 <- self$layernorm1(x + attn_output)  # (batch_size, input_seq_len, d_model)
			ffn_output <- self$ffn(out1)  # (batch_size, input_seq_len, d_model)
			ffn_output = self$dropout2(ffn_output, training = training)
			out2 <- self$layernorm2(out1 + ffn_output)  # (batch_size, input_seq_len, d_model)
			out2
		}
	})
}
