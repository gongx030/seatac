#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_autoencoder_v2_model',
		x = 'Vplots'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 100L
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		for (epoch in seq_len(epochs)){

			total_loss <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- x[b]$counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(    # convert into a C-style array
						length(b),
						x@n_intervals,
						x@n_bins_per_window,
						1L
					)) %>%
					tf$cast(tf$float32) %>%
					tf$nn$conv2d(model@gaussian_kernel, strides = c(1, 1, 1, 1), padding = 'SAME')

				xi_min <- tf$reduce_min(xi, c(1L, 2L), keepdims = TRUE)
				xi_max <- tf$reduce_max(xi, c(1L, 2L), keepdims = TRUE)
				xi <- (xi - xi_min) / (xi_max - xi_min)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					z <- xi %>% 
						model@encoder() 

					loss <- (xi - xi_pred)^2 %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()
				})

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

			flog.info(sprintf('training %s | epoch=%4.d/%4.d | total_loss=%13.7f', class(model), epoch, epochs, total_loss))

      if (epoch == 1 || epoch %% 20 == 0){
				x <- model %>% predict(x, batch_size = batch_size)
				y_umap <- umap(x$latent)$layout
				plot(y_umap, pch = 21, bg = classes, main = epoch)
			}

		}
		model
	}
)

