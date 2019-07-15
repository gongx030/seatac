fit.vae <- function(model, gr, epochs = 1, steps_per_epoch = 10, batch_size = 256, epochs_warmup = 50, beta = 1){

	optimizer <- tf$train$AdamOptimizer(0.001)
	flog.info('optimizer: Adam(learning_rate=0.001)')
	flog.info(sprintf('steps_per_epoch: %d', steps_per_epoch))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))
  num_samples <- metadata(gr)$num_samples

	batch_size <- ceiling(batch_size / num_samples)

	# determining the weight for each window
	W <- (mcols(gr)$num_reads / 50)^(3/4)
	W[W > 1] <- 1

  for (epoch in seq_len(epochs)) {

		# randomly sampling windows
		G <- matrix(sample.int(length(gr), steps_per_epoch * batch_size), steps_per_epoch, batch_size)

    total_loss <- 0
    total_loss_nll <- 0
    total_loss_kl <- 0
		beta_epoch <- min(epoch / epochs_warmup, 1) * beta

		for (s in 1:steps_per_epoch){

			b <- G[s, ]	# window index for current batch
			x <- mcols(gr)$counts[b, ] %>%
				as.matrix() %>% 
				tf$cast(tf$float32) %>% 
				tf$reshape(shape(batch_size, num_samples, model$input_dim, model$feature_dim)) %>%
				tf$reshape(shape(batch_size * num_samples, model$input_dim, model$feature_dim)) %>%
				tf$expand_dims(axis = 3L)

			if (model$batch_effect){
				g <- rep(seq_len(num_samples) - 1, batch_size) %>% 	# group index of current batch
					tf$cast(tf$int32)
			}

			w <- c(t(W[b, , drop = FALSE]))

      with(tf$GradientTape(persistent = TRUE) %as% tape, {

        approx_posterior <- x %>% model$encoder()
        approx_posterior_sample <- approx_posterior$sample()

				if (model$batch_effect){
        	decoder_likelihood <- list(approx_posterior_sample, g) %>% model$decoder()
				}else{
        	decoder_likelihood <- approx_posterior_sample %>% model$decoder()
				}

        nll <- -w * decoder_likelihood$log_prob(x)	
        avg_nll <- tf$reduce_mean(nll)

				if (model$trainable_prior)
					model$latent_prior <- model$latent_prior_model(NULL)

				if (model$prior == 'gmm'){

					kl_div <- w * (approx_posterior$log_prob(approx_posterior_sample) - model$latent_prior$log_prob(approx_posterior_sample))

				}else if (model$prior == 'hmm'){

					h <- approx_posterior_sample %>% 
						tf$reshape(shape(batch_size, num_samples, 1, latent_dim, 1)) %>% 
						tf$transpose(c(1L, 0L, 2L, 3L, 4L))

					pr <- model$latent_prior$log_prob(h) %>%
						tf$transpose(c(1L, 0L, 2L)) %>%
						tf$reshape(shape(batch_size * num_samples))
						
					kl_div <- w * (approx_posterior$log_prob(approx_posterior_sample) - pr)
				}

				kl_div <- beta_epoch * tf$reduce_mean(kl_div)
        loss <- kl_div + avg_nll
      })

      total_loss <- total_loss + loss
      total_loss_nll <- total_loss_nll + avg_nll
      total_loss_kl <- total_loss_kl + kl_div

      encoder_gradients <- tape$gradient(loss, model$encoder$variables)
      decoder_gradients <- tape$gradient(loss, model$decoder$variables)

      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients, model$encoder$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients, model$decoder$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

			if (model$trainable_prior){
      	prior_gradients <- tape$gradient(loss, model$latent_prior_model$variables)
				optimizer$apply_gradients(
					purrr::transpose(list(prior_gradients, model$latent_prior_model$variables)),
					global_step = tf$train$get_or_create_global_step()
				)
			}
    }

    flog.info(sprintf('epoch=%4.d/%4.d | negative likelihood=%7.1f | kl=%7.1f | beta=%5.3f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, beta_epoch, total_loss))
  }

} # fit_vae
