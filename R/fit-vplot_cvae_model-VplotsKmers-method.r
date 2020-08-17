#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_cvae_model',
		x = 'VplotsKmers'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 32L,
		batch_size_blocks = 256L,
		steps_per_epoch = 16L,
		epochs = 100L,
		min_reads_per_block = 10
	){

		optimizer <- tf$keras$optimizers$Adam(learning_rate)

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		for (epoch in seq_len(epochs)){

			for (j in seq_len(n_batch)){

				total_loss <- 0
				total_loss_reconstruction <- 0
				total_loss_kl <- 0

				h <- starts[j]:ends[j]

				xs <- x[h] %>% prepare_blocks(model, min_reads_per_block) 
				n_blocks <- xs$vplot$shape[0]

				for (i in seq_len(steps_per_epoch)){

					b <- sample(0:(n_blocks - 1), batch_size_blocks, replace = TRUE)
					xi <- tf$gather(xs$vplot, b)
					ci <- tf$gather(xs$kmers, b)

					with(tf$GradientTape(persistent = TRUE) %as% tape, {

						enc <- xi %>% model@encoder(ci)
						z <- enc$posterior$sample()

						xi_pred <- z %>%  
							model@decoder(enc$context)

						loss_kl <- (enc$posterior$log_prob(z) - enc$prior$log_prob(z)) %>%
							tf$reduce_mean()

						loss_reconstruction <- -xi_pred$log_prob(xi) %>%
							tf$reduce_mean()

						loss <- loss_reconstruction + loss_kl

					})

					total_loss_reconstruction  <- total_loss_reconstruction + loss_reconstruction
					total_loss_kl <- total_loss_kl + loss_kl
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

				flog.info(sprintf('training %s | epoch=%4.d/%4.d | window batch=%5.d/%5.d | n_blocks=%7.d | total_loss_reconstruction=%13.3f | total_loss_kl=%13.3f | total_loss=%13.3f', class(model), epoch, epochs, j, n_batch, n_blocks, total_loss_reconstruction, total_loss_kl, total_loss))
			}

			# evaluating the predicted performance
			valid <- sample(which(rowMeans(x$mnase) > 0.5 & rowSums(x$counts) < 25), 50)
			x_pred <- model %>% predict(x[valid])
			x_pred <- add_nucleosome_signal(x_pred)
			rmse <- sqrt(rowSums((x_pred$mnase_scaled - x_pred$nucleosome_signal)^2))
			flog.info(sprintf('evaluating %s | epoch=%4.d/%4.d | rmse=%.3f',  class(model), epoch, epochs, mean(rmse)))

		}
		model
	}
)

