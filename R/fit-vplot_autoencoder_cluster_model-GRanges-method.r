#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_autoencoder_cluster_model',
		x = 'GRanges'
	),
	function(
		model,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 100L,
		steps_per_epoch = 10L
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(model@data), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(model@data)] <- length(model@data)
		n_batch <- length(starts)

		z <- model %>% encode(batch_size = batch_size)
		km <- kmeans(z, model@num_clusters)
		u <- km$centers %>%
			tf$cast(tf$float32) %>%
			tf$Variable()

		p <- sparseMatrix(i = 1:length(model@data), j = km$cluster, dims = c(length(model@data), model@num_clusters))
		p <- p + 0.1
		p <- Diagonal(x = 1 / rowSums(p)) %*% p
		p <- p %>% 
			as.matrix() %>% 
			tf$cast(tf$float32)

		for (epoch in seq_len(epochs)){

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_cluster <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- model@data[b]$counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						length(b),
						metadata(x)$n_bins_per_window, 
						metadata(x)$n_intervals, 
						1L
					)) %>%
					tf$cast(tf$float32)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					zi <- xi %>% 
						model@encoder() 
					
					xi_pred <- zi %>%
						model@decoder()

					Pi <- p[b, ]

					di <- zi^2 %>% tf$reduce_sum(1L, keepdims = TRUE) - 
						2 * tf$matmul(zi, u, transpose_b = TRUE) + 
						tf$transpose(u)^2 %>% tf$reduce_sum(0L, keepdims = TRUE)

					loss_reconstruction <- (xi - xi_pred)^2 %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()

					loss_cluster <- (Pi * (di + model@sigma * log(Pi + 1e-3))) %>%
						tf$reduce_sum(axis = 1L) %>%
						tf$reduce_mean()
					loss_cluster <- model@gamma * loss_cluster 

					loss <- loss_reconstruction + loss_cluster

				})

				total_loss <- total_loss + loss
				total_loss_reconstruction <- total_loss_reconstruction + loss_reconstruction
				total_loss_cluster <- total_loss_cluster + loss_cluster

				encoder_gradients <- tape$gradient(loss, model@encoder$trainable_variables)
				list(encoder_gradients, model@encoder$trainable_variables) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

				decoder_gradients <- tape$gradient(loss, model@decoder$trainable_variables)
				list(decoder_gradients, model@decoder$trainable_variables) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

				other_gradients <- tape$gradient(loss, list(u))
				list(other_gradients, list(u)) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

			}

			if (epoch %% steps_per_epoch == 0){

				flog.info('updating membership')

				z <- model %>% 
					encode(batch_size = batch_size) %>%
					tf$cast(tf$float32)

				d <- z^2 %>% tf$reduce_sum(1L, keepdims = TRUE) - 
					2 * tf$matmul(z, u, transpose_b = TRUE) + 
					tf$transpose(u)^2 %>% tf$reduce_sum(0L, keepdims = TRUE)

				p <- tf$nn$softmax(-d / model@sigma)

				cls <- p %>% 
					tf$nn$softmax() %>%
					tf$math$argmax(axis = 1L) %>%
					as.numeric()
				cls <- cls + 1

				z %>%
					as.matrix() %>%
					umap() %>%
					pluck('layout') %>%
					plot(pch = 21, bg = cls, col = cls, main = epoch, cex = 0.25)

#				y <- u %>% 
#					model@decoder()

				for (k in 1:model@num_clusters){

					if (sum(cls == k) > 1){
						model@data$counts[cls == k, ] %>% 
							colMeans() %>% 
							matrix(metadata(model@data)$n_bins_per_window , metadata(model@data)$n_intervals) %>% 
							image(main = k)
					}
				}
				print(table(cls, model@data$sample_id))
			}

			flog.info(sprintf('training vplot_autoencoder_cluster_model | epoch=%4.d/%4.d | loss_reconstruction=%13.7f | loss_cluster=%13.7f | total_loss=%13.7f', epoch, epochs, total_loss_reconstruction, total_loss_cluster, total_loss))

		}

		model@membership <- as.matrix(p)

		model@centers <- u %>%
			tf$convert_to_tensor() %>%
			as.matrix()

		model

	}
)

