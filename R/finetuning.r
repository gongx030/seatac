#' finetuning
#'
finetuning <- function(model, centers = 2, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	window_dim <- length(gr)  # number of samples

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
  flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

	for (epoch in seq_len(epochs)){

		# k-means
		gr <- model %>% predict(model$data)
		km <- kmeans(gr$latent, centers)
		mu <- sparseMatrix(i = 1:window_dim, km$cluster, dims = c(window_dim, centers)) %*% km$centers %>%
			as.matrix()

		total_loss <- 0

		# updating networks 
		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			fragment_size <- mcols(gr[b])$fragment_size %>%
				as.matrix() %>%
				tf$cast(tf$float32)

			position <- mcols(gr[b])$position %>%
				as.matrix() %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- list(fragment_size, position) %>% model$encoder()
				posterior_sample <- posterior$sample()

				loss <- tf$reduce_mean(tf$reduce_sum(tf$square(posterior_sample, mu), 1L))
			})

			total_loss <- total_loss + loss

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
		}

		flog.info(sprintf('finetuning | epoch=%4.d/%4.d | total=%7.1f', epoch, epochs, total_loss))

		gr <- model %>% predict(model$data)
		y <- Rtsne(gr$latent, check_duplicates = FALSE)$Y
		smoothScatter(y, main = 'v-plot')
#		km$centers %>% points(pch = 21, bg = 'red')
	}

	browser()
} # finetuning
