#' fit
#'
fit.vae_baseline <- function(model, ...) fit_vae_vplot(model, ...)


#' fit_vae_vplot
#'
fit_vae_vplot <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

	if (length(model$prior$trainable_variables) > 0)
		trainable_prior <- TRUE
	else
		trainable_prior <- FALSE

	flog.info(sprintf('trainable prior: %s', trainable_prior))

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x)
				avg_nll <- tf$reduce_mean(nll)

				kl_div <- posterior$log_prob(posterior_sample) - model$prior(NULL)$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 
			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			if (trainable_prior){
				prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)
			}

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

			if (trainable_prior){
				optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))
			}

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit_vae_vplot
