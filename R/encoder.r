#' encoder_model
#'
#' @param window_dim window dimension of the v-plot (i.e. the width of the genomic region)
#' @param interval_dim interval dimension of the v-plot	(i.e. the fragment size dimension)
#' @param filters0 the beginning filter dimension coming out of the latent layer
#' @param filters the filter sizes of each deconv layer
#' @param kernel_size the kernel size of each deconv layer.  The interval and window spaces shared the same kernel size. 
#'
encoder_model <- function(
	sample_embedding_dim = 20L,
	sequence_embedding_dim = 100L,
	sample_dim,
	latent_dim = 50L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$embedding_1 <- layer_embedding(
			input_dim = 4L,
			output_dim = 4L,
			weights = list(diag(4L)),
			trainable = FALSE
		) 

		self$sequence_conv <- layer_conv_1d(
			filters = sequence_embedding_dim,
			kernel_size = 6L,
			strides = 1L,
			activation = 'relu'
		)

		self$bn_1 <- layer_batch_normalization()

		self$bigru <- bidirectional(layer = layer_cudnn_gru(units = 50L))

		self$embedding_2 <- layer_embedding(
			input_dim = sample_dim,
			output_dim = sample_embedding_dim
		)

		self$dropout_1 <- layer_dropout(rate = 0.2)

		self$dense_1 <- layer_dense(
			units = latent_dim,
			activation = 'relu'
		)

		function(x, mask = NULL){

			y_seq <- x[[1]] %>%
				self$embedding_1() %>%
				self$sequence_conv() %>%
				self$bn_1() %>%
				self$bigru()

			y_sample <- x[[2]] %>%
				self$embedding_2()

			y <- list(y_seq, y_sample) %>%
				layer_concatenate() %>%
				self$dropout_1() %>%
				self$dense_1()

			y

		}
	})
} # decoder_model_vae_baseline_conv


