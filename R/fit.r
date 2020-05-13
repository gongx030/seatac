
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
		recurrent_steps = 3,
		K = 10L	# number of samples drawn from each V-plot
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))
		flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		for (epoch in seq_len(epochs)) {

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_kl <- 0

			for (s in 1:steps_per_epoch){

				b <- sample.int(length(x), batch_size)

				xi <- mcols(x[b])$counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						batch_size, 
						metadata(x)$n_bins_per_window, metadata(x)$n_intervals, 
						1L
					)) %>%
					tf$cast(tf$float32)

				w <- (xi > 0) %>%	# index for the non-zero terms
					tf$cast(tf$float32)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					for (iter in seq_len(recurrent_steps)){

						if (iter == 1)
							xi_input <- xi
						else
							xi_input <- xi * w + xi_pred * (1 - w)

						posterior <- xi_input %>% model@encoder()
						posterior_sample <- posterior$sample()

						xi_pred <- posterior_sample %>% model@decoder()

						xi_pred <- xi_pred$mean() %>%
							tf$clip_by_value(clip_value_min = 0, clip_value_max = 100)

					}

					print(range(as.matrix(posterior$mean())))

					print(range(as.matrix(xi_pred)))

#					loss_reconstruction <- (w * (xi - xi_pred)^2) %>%
					loss_reconstruction <- ((xi - xi_pred)^2) %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()

					kl_div <- (posterior$log_prob(posterior_sample) - model@prior(NULL)$log_prob(posterior_sample)) %>%
						tf$reduce_mean()

					loss <- loss_reconstruction  + kl_div

				})

				total_loss <- total_loss + loss
				total_loss_reconstruction <- total_loss_reconstruction + loss_reconstruction
				total_loss_kl <- total_loss_kl + kl_div

				encoder_gradients <- tape$gradient(loss, model@encoder$trainable_variables)
				list(encoder_gradients, model@encoder$trainable_variables) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

				decoder_gradients <- tape$gradient(loss, model@decoder$trainable_variables)
				list(decoder_gradients, model@decoder$trainable_variables) %>%
					purrr::transpose() %>%
					optimizer$apply_gradients()

			}


			flog.info(sprintf('training | epoch=%4.d/%4.d | loss_reconstruction=%9.1f | kl=%9.1f | total_loss=%9.1f', epoch, epochs, total_loss_reconstruction, total_loss_kl, total_loss))

			if (epoch == 1 || epoch %% 5 == 0){

				ii <- sample(1:length(x), 2000)

				z <- model %>% encode(x[ii], batch_size = batch_size, recurrent_steps = recurrent_steps)
				X <- x$counts[ii, ]
				sample_label <- x$sample_id[ii]

				kk <- 3

				y_umap <- umap(z)$layout
				cls <- kmeans(z, kk, nstart = 10)$cluster
				y_umap %>% plot(pch = 21, bg = sample_label)
				y_umap %>% plot(pch = 21, bg = cls)

				for (jj in 1:kk){
					X[cls == jj, ] %>% 
						colMeans() %>% 
						matrix(metadata(x)$n_bins_per_window , metadata(x)$n_intervals) %>% 
						image(main = jj)
				}
				print(table(cls, sample_label))
			}
		}
	}
) # fit
