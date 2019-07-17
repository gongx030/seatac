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
    total_loss_nll_x <- 0
    total_loss_nll_y <- 0
    total_loss_kl_x <- 0
    total_loss_kl_y <- 0

		for (s in 1:steps_per_epoch){

			b <- G[s, ]	# window index for current batch
			x <- mcols(gr)$counts[b, ] %>%
				as.matrix() %>% 
				tf$cast(tf$float32) %>% 
				tf$reshape(shape(batch_size, num_samples, model$input_dim, model$feature_dim)) %>%
				tf$reshape(shape(batch_size * num_samples, model$input_dim, model$feature_dim)) %>%
				tf$expand_dims(axis = 3L)
			
			y <- mcols(gr)$coverage[b, , , drop = FALSE] %>%
				tf$cast(tf$float32) %>% 
				tf$transpose(c(0L, 2L, 1L)) %>% 
				tf$reshape(shape(batch_size * num_samples, model$feature_dim)) %>%
				tf$expand_dims(axis = 2L)

			w <- c(t(W[b, , drop = FALSE]))

      with(tf$GradientTape(persistent = TRUE) %as% tape, {

        posterior_x <- x %>% model$encoder$vplot()
        posterior_sample_x <- posterior_x$sample()

        posterior_y <- y %>% model$encoder$coverage()
        posterior_sample_y <- posterior_y$sample()

				posterior_sample <- tf$concat(list(posterior_sample_x, posterior_sample_y), 1L)

       	likelihood_x <- posterior_sample %>% model$decoder$vplot()
       	likelihood_y <- posterior_sample %>% model$decoder$coverage()

				nll_x <- -w * likelihood_x$log_prob(x)
				nll_y <- -w * likelihood_y$log_prob(y)

        avg_nll_x <- tf$reduce_mean(nll_x)
        avg_nll_y <- tf$reduce_mean(nll_y)

				vplot_prior <- model$latent_prior_model$vplot(NULL)

				if (model$prior == 'gmm'){

					kl_div_x <- w * (posterior_x$log_prob(posterior_sample_x) - vplot_prior$log_prob(posterior_sample_x))

				}else if (model$prior == 'hmm'){

#					h <- approx_posterior_sample %>% 
#						tf$reshape(shape(batch_size, num_samples, 1, latent_dim, 1)) %>% 
#						tf$transpose(c(1L, 0L, 2L, 3L, 4L))

#					pr <- model$latent_prior$log_prob(h) %>%
#						tf$transpose(c(1L, 0L, 2L)) %>%
#						tf$reshape(shape(batch_size * num_samples))
						
#					kl_div <- w * (approx_posterior$log_prob(approx_posterior_sample) - pr)
				}

				kl_div_x <- beta * tf$reduce_mean(kl_div_x)

				coverage_prior <- model$latent_prior_model$coverage(NULL)
				kl_div_y <- w * (posterior_y$log_prob(posterior_sample_y) - coverage_prior$log_prob(posterior_sample_y))
				kl_div_y <- beta * tf$reduce_mean(kl_div_y)

        loss <- kl_div_x + kl_div_y + avg_nll_x + avg_nll_y
      })

      total_loss <- total_loss + loss
			total_loss_nll_x <- total_loss_nll_x + avg_nll_x
      total_loss_nll_y <- total_loss_nll_y + avg_nll_y
			total_loss_kl_x <- total_loss_kl_x + kl_div_x
      total_loss_kl_y <- total_loss_kl_y + kl_div_y

      encoder_gradients_x <- tape$gradient(loss, model$encoder$vplot$variables)
      encoder_gradients_y <- tape$gradient(loss, model$encoder$coverage$variables)
      decoder_gradients_x <- tape$gradient(loss, model$decoder$vplot$variables)
      decoder_gradients_y <- tape$gradient(loss, model$decoder$coverage$variables)

      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients_x, model$encoder$vplot$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients_y, model$encoder$coverage$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients_x, model$decoder$vplot$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients_y, model$decoder$coverage$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

     	prior_gradients <- tape$gradient(loss, model$latent_prior_model$vplot$variables)
			optimizer$apply_gradients(
				purrr::transpose(list(prior_gradients, model$latent_prior_model$vplot$variables)),
				global_step = tf$train$get_or_create_global_step()
			)
    }

    flog.info(sprintf('epoch=%4.d/%4.d | nll(vplot)=%7.1f | nll(coverage)=%7.1f | kl(vplot)=%7.1f | kl(coverage)=%7.1f | beta=%5.3f | total=%7.1f', epoch, epochs, total_loss_nll_x, total_loss_nll_y, total_loss_kl_x, total_loss_kl_y, beta, total_loss))
  }

} # fit_vae
