#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_vae_model',
		x = 'Vplots'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 32L,
		batch_size_blocks = 256L,
		steps_per_epoch = 16L,
		epochs = 100L
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

				xs <- x[h] %>% prepare_blocks(model)

				for (i in seq_len(steps_per_epoch)){

					b <- sample(0:(xs$shape[0] - 1), batch_size_blocks, replace = TRUE)
					xi <- tf$gather(xs, b)

					with(tf$GradientTape(persistent = TRUE) %as% tape, {

						posterior <- xi %>%
							model@encoder()

						z <- posterior$sample()

						xi_pred <- z %>%  
							model@decoder()

						loss_reconstruction <- -xi_pred$log_prob(xi) %>%
							tf$reduce_mean()

						loss_kl <- (posterior$log_prob(z) - model@prior(NULL)$log_prob(z)) %>%
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

				flog.info(sprintf('training %s | epoch=%4.d/%4.d | window batch=%5.d/%5.d | total_loss_reconstruction=%13.3f | total_loss_kl=%13.3f | total_loss=%13.3f', class(model), epoch, epochs, j, n_batch, total_loss_reconstruction, total_loss_kl, total_loss))
			}

		}
		model
	}
)

