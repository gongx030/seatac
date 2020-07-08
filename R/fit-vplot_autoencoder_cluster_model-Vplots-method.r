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
		x,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 100L,
		epochs_per_updating = 10L
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		# initialize the clusters
		X <- x$counts %>%
		  as.matrix() %>%
		    array(c(
					length(x),
					x@n_bins_per_window,
					x@n_intervals
			))

		X_window <- rowSums(X, dims = 2)
		X_window <- Diagonal(x = 1 / rowSums(X_window)) %*% X_window

		X <- aperm(X, c(1, 3, 2))
		X_interval <- rowSums(X, dims = 2)
		X_interval <- Diagonal(x = 1 / rowSums(X_interval)) %*% X_interval

		km <- kmeans(cbind(X_window, X_interval), nstart = 10, centers = model@num_clusters)
		print(table(km$cluster, classes))

		# initialize the centers
		u <- tf$Variable(tf$random$normal(shape(model@num_clusters, model@latent_dim)))

		# memership matrix
		p <- sparseMatrix(i = 1:length(x), j = km$cluster, dims = c(length(x), model@num_clusters))
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

				xi <- x[b]$counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
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

			if (epoch %% epochs_per_updating == 0){

				flog.info('updating membership')

				x <- model %>% predict(x, batch_size = batch_size)

				z <- x$latent %>% 
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

				
				print(table(cls, classes))

			}

			if (epoch == 1 || epoch %% 10 == 0){
				x <- model %>% predict(x, batch_size = batch_size)
				y_umap <- umap(x$latent)$layout
#				print(table(kmeans(y_umap, 3)$cluster, classes))
				plot(y_umap, pch = 21, bg = classes, col = classes, main = epoch, cex = 0.5)
			}

			flog.info(sprintf('training %s | epoch=%4.d/%4.d | loss_reconstruction=%13.7f | loss_cluster=%13.7f | total_loss=%13.7f', class(model), epoch, epochs, total_loss_reconstruction, total_loss_cluster, total_loss))

		}

		model@membership <- as.matrix(p)

		model@centers <- u %>%
			tf$convert_to_tensor() %>%
			as.matrix()

		model

	}
)

