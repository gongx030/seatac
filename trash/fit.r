#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_autoencoder_model2',
		x = 'GRanges'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 50L,
		burnin = 10L
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))
		flog.info(sprintf('beta(beta):%.3e', beta))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		for (epoch in seq_len(burnin)){

			total_loss <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- x[b]$counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						length(b),
						metadata(x)$n_bins_per_window, 
						metadata(x)$n_intervals, 
						1L
					)) %>%
					tf$cast(tf$float32)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					xi_pred <- xi %>% 
						model@encoder() %>%
						model@decoder()

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

			flog.info(sprintf('burnin | epoch=%4.d/%4.d | total_loss=%9.1f', epoch, burnin, total_loss))

		}

		browser()

		u <- tf$random$normal(shape(model@num_clusters, model@latent_dim)) %>%
			tf$Variable(dtype = tf$float32)

		g <- tf$random$normal(shape(length(x), model@num_clusters)) %>%
			tf$Variable(dtype = tf$float32)

		for (epoch in seq_len(epochs)){

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_cluster <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- x[b]$counts  %>%
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

					gi <- g[b, ] %>% 
						tf$nn$softmax()

#					di <- xi^2 %>% tf$reduce_sum(c(1L, 2L)) - 
#						2 * tf$matmul(
#							xi %>% tf$reshape(shape(length(b), metadata(x)$n_bins_per_window * metadata(x)$n_intervals)),
#							yi_pred %>% tf$reshape(shape(model@num_clusters, metadata(x)$n_bins_per_window * metadata(x)$n_intervals)),
#							transpose_b = TRUE
#						)	+
#						yi_pred^2 %>% 
#							tf$reduce_sum(c(1L, 2L)) %>%
#							tf$transpose()

					di <- zi^2 %>% tf$reduce_sum(1L, keepdims = TRUE) - 
						2 * tf$matmul(zi, u, transpose_b = TRUE) + 
						tf$transpose(u)^2 %>% tf$reduce_sum(0L, keepdims = TRUE)

					loss_reconstruction <- (xi - xi_pred)^2 %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()

					loss_cluster <- (gi * (di + log(gi + 1e-3))) %>%
						tf$reduce_sum(axis = 1L) %>%
						tf$reduce_mean()
					loss_cluster <- beta * loss_cluster 

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

				others_gradients <- tape$gradient(loss, list(g, u))
				list(others_gradients, list(g, u)) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

			}

			flog.info(sprintf('training | epoch=%4.d/%4.d | loss_reconstruction=%9.1f | loss_cluster=%9.1f | total_loss=%9.1f', epoch, epochs, total_loss_reconstruction, total_loss_cluster, total_loss))

			if (epoch %% 10 == 0){

				ii <- sample(1:length(x), 1000)

				zii <- model %>% 
					encode(x[ii], batch_size = batch_size)

				cls <- g %>% 
					tf$nn$softmax() %>%
					tf$math$argmax(axis = 1L) %>%
					as.numeric()
				cls <- cls[ii]
				cls <- cls + 1

				Xii <- x$counts[ii, ]
				sample_label <- x$sample_id[ii]

#				yii <- u %>%
#					model@decoder()

				y_umap <- umap(zii)$layout

				y_umap %>% plot(pch = 21, bg = sample_label)
				y_umap %>% plot(pch = 21, bg = cls)

				for (jj in 1:model@num_clusters){
					Xii[cls == jj, ] %>% 
						colMeans() %>% 
						matrix(metadata(x)$n_bins_per_window , metadata(x)$n_intervals) %>% 
						image(main = jj)
#					yii[jj, , , 1] %>%
#						as.matrix() %>%
#						t() %>%
#						image(main = jj)
				}
				print(table(cls, sample_label))
			}
		}
	}
) # fit


#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_knn_autoencoder_model',
		x = 'GRanges'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 50L,
		burnin = 10L,
		beta = 1
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))
		flog.info(sprintf('beta(beta):%.3e', beta))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		for (epoch in seq_len(burnin)){

			total_loss <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- x[b]$counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						length(b),
						metadata(x)$n_bins_per_window, 
						metadata(x)$n_intervals, 
						1L
					)) %>%
					tf$cast(tf$float32)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					xi_pred <- xi %>% 
						model@encoder() %>%
						model@decoder()

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

			flog.info(sprintf('burnin | epoch=%4.d/%4.d | total_loss=%9.1f', epoch, burnin, total_loss))

		}

		z <- model %>% encode(x, batch_size = batch_size)
		km <- kmeans(z, model@num_clusters, nstart = 10)

		u <- km$centers %>%
			tf$Variable(dtype = tf$float32)

		g <- sparseMatrix(i = 1:length(x), j = km$cluster, dims = c(length(x), model@num_clusters))
		g <- log(g + 0.1) %>%
			as.matrix() %>%
			tf$Variable(dtype = tf$float32)

		for (epoch in seq_len(epochs)){

			z <- model %>% encode(x, batch_size = batch_size)
			G <- FNN::knn.index(z, model@K)
			G <- sparseMatrix(
				i = rep(1:length(x), model@K), 
				j = c(G), 
				x = rep(1 / model@K, length(x) * model@K),
				dims = c(length(x), length(x))
			)

			X <- G %*% x$counts

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_cluster <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- X[b, ] %>%
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

					gi <- g[b, ] %>% 
						tf$nn$softmax()

					di <- zi^2 %>% tf$reduce_sum(1L, keepdims = TRUE) - 
						2 * tf$matmul(zi, u, transpose_b = TRUE) + 
						tf$transpose(u)^2 %>% tf$reduce_sum(0L, keepdims = TRUE)

					xi_pred <- zi %>%
						model@decoder()

					loss_reconstruction <- (xi - xi_pred)^2 %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()

					loss_cluster <- (gi * (di + log(gi + 1e-3))) %>%
						tf$reduce_sum(axis = 1L) %>%
						tf$reduce_mean()
					loss_cluster <- beta * loss_cluster 

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

				others_gradients <- tape$gradient(loss, list(g, u))
				list(others_gradients, list(g, u)) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

			}

			flog.info(sprintf('training | epoch=%4.d/%4.d | loss_reconstruction=%9.1f | loss_cluster=%9.1f | total_loss=%9.1f', epoch, epochs, total_loss_reconstruction, total_loss_cluster, total_loss))

			if (epoch == 1 || epoch %% 10 == 0){

				ii <- sample(1:length(x), 1000)
				Xii <- X[ii, ]
				sample_label <- x$sample_id[ii]
				zii <- model %>% encode(x[ii], batch_size = batch_size)
				cls <- g %>% 
					tf$nn$softmax() %>%
					tf$math$argmax(axis = 1L) %>%
					as.numeric()
				cls <- cls[ii]
				cls <- cls + 1

				zii %>% plot(pch = 21, bg = cls, col = cls)
				u %>% tf$convert_to_tensor() %>% as.matrix() %>% points(pch = 2, col = cls, cex = 3)

				yii <- u %>%
					model@decoder()

				for (jj in 1:model@num_clusters){
					if (sum(cls == jj) > 1){
						Xii[cls == jj, ] %>% 
							colMeans() %>% 
							matrix(metadata(x)$n_bins_per_window , metadata(x)$n_intervals) %>% 
							image(main = jj)
						yii[jj, , , 1] %>%
							as.matrix() %>%
							t() %>%
							image(main = jj)
					}
				}
				print(table(cls, sample_label))
			}

		}
	}
)


