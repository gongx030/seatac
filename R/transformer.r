TransformerEncoderLayer <- reticulate::PyClass(
	'TransformerEncoderLayer',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(self, d_model, num_heads, dff, rate = 0.1) {

			super()$`__init__`()

			self$mha <- MultiHeadAttention(d_model, num_heads)
			self$ffn <- point_wise_feed_forward_network(d_model, dff)

			self$layernorm1 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
			self$layernorm2 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
						    
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


TransformerDecoderLayer <- reticulate::PyClass(
	'TransformerDecoderLayer',
	inherit = tf$keras$layers$Layer,
	list(
		`__init__` = function(self, d_model, num_heads, dff, rate = 0.1) {

			super()$`__init__`()

			self$mha1 <- MultiHeadAttention(d_model, num_heads)
			self$mha2 <- MultiHeadAttention(d_model, num_heads)

			self$ffn <- point_wise_feed_forward_network(d_model, dff)

			self$layernorm1 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
			self$layernorm2 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
			self$layernorm3 <- tf$keras$layers$LayerNormalization(epsilon = 1e-6)
						    
			self$dropout1 <- tf$keras$layers$Dropout(rate)
			self$dropout2 <- tf$keras$layers$Dropout(rate)
			self$dropout3 <- tf$keras$layers$Dropout(rate)

			NULL

		},
		call = function(self, x, enc_output, training, look_ahead_mask, padding_mask){

			res <- self$mha1(x, x, x, look_ahead_mask)  # (batch_size, target_seq_len, d_model)
			attn1 <- res$output
			attn_weights_block1 <- res$attention_weights

			attn1 <- self$dropout1(attn1, training = training)
		  out1 <- self$layernorm1(attn1 + x)
				    
			res <- self$mha2(enc_output, enc_output, out1, padding_mask)  # (batch_size, target_seq_len, d_model)
	    attn2 <- res$output
			attn_weights_block2 <- res$attention_weights
																									         
			attn2 <- self$dropout2(attn2, training = training)
			out2 <- self$layernorm2(attn2 + out1)  # (batch_size, target_seq_len, d_model)
						    
			ffn_output <- self$ffn(out2)  # (batch_size, target_seq_len, d_model)
			ffn_output <- self$dropout3(ffn_output, training = training)
			out3 <- self$layernorm3(ffn_output + out2)  # (batch_size, target_seq_len, d_model)
										    
			list(
				output = out3, 
				attention_weights_block1 = attn_weights_block1, 
				attention_weights_block2 = attn_weights_block2
			)
		}
	)
)
