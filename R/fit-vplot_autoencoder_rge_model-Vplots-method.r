#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_autoencoder_rge_model'
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

		w <- diag(model@num_clusters) %>%
			tf$cast(tf$float32)

		for (epoch in seq_len(epochs)){

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_cluster <- 0
			total_loss_graph_embedding <- 0

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

					v <- u %>% 
						model@decoder() %>%
						tf$reshape(shape(model@num_clusters, metadata(model@data)$n_bins_per_window * metadata(model@data)$n_intervals))

					v2 <- v^2 %>% tf$reduce_sum(1L, keepdims = TRUE)
					dv <- v2 - 2 * tf$matmul(v, v, transpose_b = TRUE) + tf$transpose(v2)

					loss_reconstruction <- (xi - xi_pred)^2 %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()

					loss_cluster <- (Pi * (di + model@sigma * log(Pi + 1e-3))) %>%
						tf$reduce_sum(axis = 1L) %>%
						tf$reduce_mean()
					loss_cluster <- model@gamma * loss_cluster 

					loss_graph_embedding <- (w * dv) %>%
						tf$reduce_sum(axis = shape(1L)) %>%
						tf$reduce_mean()
					loss_graph_embedding <- model@lambda * loss_graph_embedding

					loss <- loss_reconstruction + loss_cluster + loss_graph_embedding

				})

				total_loss <- total_loss + loss
				total_loss_reconstruction <- total_loss_reconstruction + loss_reconstruction
				total_loss_cluster <- total_loss_cluster + loss_cluster
				total_loss_graph_embedding <- total_loss_graph_embedding + loss_graph_embedding

				encoder_gradients <- tape$gradient(loss, model@encoder$trainable_variables)
				list(encoder_gradients, model@encoder$trainable_variables) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

				decoder_gradients <- tape$gradient(loss, model@decoder$trainable_variables)
				list(decoder_gradients, model@decoder$trainable_variables) %>%
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

				flog.info('updating graph')

				g <- u %>%
					model@decoder() %>%
					tf$reshape(shape(model@num_clusters, metadata(model@data)$n_bins_per_window * metadata(model@data)$n_intervals)) 

				g2 <- g^2 %>% tf$reduce_sum(1L, keepdims = TRUE)
				dg <- g2 - 2 * tf$matmul(g, g, transpose_b = TRUE) + tf$transpose(g2)
				dg <- dg / (dg + 1)	# distance to similarity

				w <- dg  %>% 
					as.matrix() %>%
					as('dgCMatrix') %>%
					igraph::graph_from_adjacency_matrix(weighted = TRUE) %>%
					igraph::mst() %>%
					igraph::as_adjacency_matrix() %>%
					as.matrix() %>%
					tf$cast(tf$float32)
				w <- w + tf$transpose(w)

				flog.info('updating centers')

				l <- w %>% tf$reduce_sum(1L) %>% tf$linalg$diag() - w	# graph laplacian matrix
				L <- tf$linalg$diag(tf$reduce_sum(p, axis = 0L))

				u <- tf$matmul(
					tf$matmul(z, p, transpose_a = TRUE),
					(2 * model@lambda * l + model@gamma * L) %>%
						tf$linalg$inv()
				) %>%
					tf$transpose()
				
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

				for (k in 1:model@num_clusters){

#					if (sum(cls == k) > 1){
#						model@data$counts[cls == k, ] %>% 
#							colMeans() %>% 
#							matrix(metadata(model@data)$n_bins_per_window , metadata(model@data)$n_intervals) %>% 
#							image(main = k)
#					}
				}
				print(table(cls, model@data$sample_id))
			}

			flog.info(sprintf('training vplot_autoencoder_rge_model | epoch=%4.d/%4.d | loss_reconstruction=%13.7f | loss_cluster=%13.7f | loss_graph_embedding=%13.7f | total_loss=%13.7f', epoch, epochs, total_loss_reconstruction, total_loss_cluster, total_loss_graph_embedding, total_loss))

		}

		model@membership <- as.matrix(p)

		model@centers <- u %>%
			tf$convert_to_tensor() %>%
			as.matrix()

		model

	}
)

