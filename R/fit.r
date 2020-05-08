
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
		steps_per_epoch = 10
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))
		flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		interval_per_batch <- floor(batch_size / metadata(x)$n_samples)	
		batch_size2 <- interval_per_batch * metadata(x)$n_samples

		for (epoch in seq_len(epochs)) {

			total_loss <- 0
			total_loss_reconstruction <- 0
			total_loss_kl <- 0

			for (s in 1:steps_per_epoch){

				b <- sample.int(length(x), interval_per_batch)

				xs <- mcols(x[b])$counts %>%
					array_reshape(c(
						length(b), 
						metadata(x)$n_bins_per_window,
						metadata(x)$n_intervals,
						metadata(x)$n_samples
						)
					) %>%
					array_permute(c(1, 4, 2, 3)) %>%
					array_reshape(c(
						batch_size2,
						metadata(x)$n_bins_per_window *metadata(x)$n_intervals
					)) %>%
					as_dgCMatrix() %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						batch_size2, 
						metadata(x)$n_bins_per_window, metadata(x)$n_intervals, 
						1L
					)) %>%
					tf$cast(tf$float32)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					posterior <- xs %>% model@encoder()

					posterior_sample <- posterior$sample()

					likelihood <- posterior_sample %>% model@decoder()

					loss_reconstruction <- -likelihood$log_prob(xs) %>%
						tf$reduce_mean()

					kl_div <- (posterior$log_prob(posterior_sample) - model@prior(NULL)$log_prob(posterior_sample)) %>%
						tf$reduce_mean()

					loss <- loss_reconstruction + kl_div
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

#			xs %>% tf$reduce_sum(0L) %>% tf$squeeze() %>% as.matrix() %>% image()
#			likelihood$mean() %>% tf$reduce_sum(0L) %>% tf$squeeze() %>% as.matrix() %>% image()

			flog.info(sprintf('training | epoch=%4.d/%4.d | loss_reconstruction=%9.1f | kl=%9.1f | total_loss=%9.1f', epoch, epochs, total_loss_reconstruction, total_loss_kl, total_loss))

		}

	}
) # fit
