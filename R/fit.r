
#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_autoencoder_model',
		x = 'GRanges'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 128, 	# v-plot per batch
		epochs = 50, 
		steps_per_epoch = 10,
		K = 8L # number of samples drawn from each V-plot
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))
		flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		for (epoch in seq_len(epochs)) {

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_center_reconstruction <- 0
			total_loss_cluster <- 0

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

					zi <- xi %>% 
						model@encoder() 

					xi_pred <- zi %>%
						model@decoder()

					di <- zi^2 %>% tf$reduce_sum(1L, keepdims = TRUE) - 
						2 * tf$matmul(zi, u, transpose_b = TRUE) + 
						tf$transpose(u)^2 %>% tf$reduce_sum(0L, keepdims = TRUE)

					Pi <- tf$nn$softmax(-di)

					yi <- tf$tensordot(tf$transpose(Pi), xi, 1L) / tf$cast(batch_size, tf$float32)

					yi_pred <- u %>%
						model@decoder()

					loss_cluster <- (Pi * (di + log(Pi + 1e-3))) %>%
						tf$reduce_sum(axis = 1L) %>%
						tf$reduce_mean()
					
					loss_reconstruction <- loss_mean_squared_error(xi, xi_pred) %>%
						tf$reduce_sum(axis = shape(1L, 2L)) %>%
						tf$reduce_mean()

					loss_center_reconstruction <- loss_mean_squared_error(yi, yi_pred) %>%
						tf$reduce_sum(axis = shape(1L, 2L)) %>%
						tf$reduce_mean()

					loss <- loss_reconstruction  + loss_cluster + loss_center_reconstruction

				})

				total_loss <- total_loss + loss
				total_loss_reconstruction <- total_loss_reconstruction + loss_reconstruction
				total_loss_center_reconstruction <- total_loss_center_reconstruction + loss_center_reconstruction
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

			flog.info(sprintf('training | epoch=%4.d/%4.d | loss_reconstruction=%9.1f | loss_center_reconstruction=%9.1f | total_loss_cluster=%9.1f | total_loss=%9.1f', epoch, epochs, total_loss_reconstruction, total_loss_center_reconstruction, total_loss_cluster, total_loss))

#			P <- tf$nn$softmax(g) %>%
#				as.matrix()
#			U <- t(P) %*% x$counts

#			if (epoch == 5) browser()

			print(u)

			for (kk in 1:K){
				y <- u %>%
					model@decoder()
				y[kk, , , 1] %>% as.matrix() %>% t() %>% image(main = kk)
			}

			if (epoch > 100000){

				ii <- sample(1:length(x), 2000)

				z <- model %>% encode(x[ii], batch_size = batch_size)
				Xii <- x$counts[ii, ]
				sample_label <- x$sample_id[ii]

				kk <- 3

				y_umap <- umap(z)$layout
				cls <- kmeans(y_umap, kk, nstart = 10)$cluster
				y_umap %>% plot(pch = 21, bg = sample_label)
				y_umap %>% plot(pch = 21, bg = cls)

				for (jj in 1:kk){
					Xii[cls == jj, ] %>% 
						colMeans() %>% 
						matrix(metadata(x)$n_bins_per_window , metadata(x)$n_intervals) %>% 
						image(main = jj)
				}
				print(table(cls, sample_label))
			}
		}
	}
) # fit
