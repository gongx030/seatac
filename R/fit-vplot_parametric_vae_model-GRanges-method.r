#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_parametric_vae_model',
		x = 'GRanges'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 100L
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		n <- rowSums(x$counts)	# total number of reads per V-plot

		for (epoch in seq_len(epochs)){

			total_loss <- 0

			total_loss_likelihood <- 0
			total_loss_kl <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				# data input
				xi <- x[b]$smoothed_counts %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						length(b),
						metadata(x)$n_intervals, 
						metadata(x)$n_bins_per_window, 
						1L
					)) %>%
					tf$cast(tf$float32)

#				y <- x[b]$counts %>%
				y <- x[b]$smoothed_counts %>%
					array_reshape(c(	
						length(b),
						metadata(x)$n_bins_per_window,
						metadata(x)$n_intervals
					))

				m <- y@vals %>% 
					tf$cast(tf$float32) %>%
					tf$reshape(c(length(y@vals), 1L))

				y <- y@subs

				g <- sparseMatrix(i = 1:nrow(y), j = y[, 1], dims = c(nrow(y), length(b))) %>%
					as.matrix() %>%
					tf$cast(tf$float32)

				# read position relative to the center
				pos <- metadata(x)$positions[y[, 2]] %>%
					abs() %>%
					tf$cast(tf$float32) %>%
					tf$reshape(c(nrow(y), 1L))

				# probability of NFR
				p <- model@nfr_prob[y[, 3]] %>%
					tf$cast(tf$float32) %>%
					tf$reshape(c(nrow(y), 1L))

				# weight for each v-plot
				w <- 1 / n[y[, 1]]	%>%	
					tf$cast(tf$float32) %>%
					tf$reshape(c(nrow(y), 1L))

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					posterior <- xi %>% 
						model@encoder()

					z <- posterior$sample()

					u <- z %>% model@decoder()

					mu_nfr <- tf$matmul(g, u[, 1, drop = FALSE])

					sd_nfr <- tf$matmul(g, u[, 2, drop = FALSE])
					sd_nfr <- sd_nfr + 1

					mu_mono <- tf$matmul(g, u[, 3, drop = FALSE])

					sd_mono<- tf$matmul(g, u[, 4, drop = FALSE])
					sd_mono <- sd_mono + 1

					loss_likelihood <- (-w * m * (p * tfd_normal(mu_nfr, sd_nfr)$log_prob(pos) + (1 - p) * tfd_normal(mu_mono, sd_mono)$log_prob(pos))) %>%
						tf$reduce_mean()

					loss_kl <- (posterior$log_prob(z) - model@prior(NULL)$log_prob(z)) %>%
						tf$reduce_mean()

					loss <- loss_likelihood + loss_kl 

				})


				total_loss_likelihood  <- total_loss_likelihood + loss_likelihood
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

				print(mu_nfr %>% as.numeric() %>% range())
			flog.info(sprintf('training %s | epoch=%4.d/%4.d | total_loss_likelihood=%13.7f | total_loss_kl=%13.7f | total_loss=%13.7f', class(model), epoch, epochs, total_loss_likelihood, total_loss_kl, total_loss))

		}


		x0 <- x$smoothed_counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(		# convert into a C-style array
				length(x),
				metadata(x)$n_intervals, 
				metadata(x)$n_bins_per_window, 
				1L
			)) %>%
			tf$cast(tf$float32)

#		x$smoothed_counts %>% as.matrix() %>% colMeans() %>% matrix(metadata(x)$n_bins_per_window, metadata(x)$n_intervals) %>% image()
#		x0 %>% tf$reduce_mean(0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image()
#		xi %>% tf$reduce_mean(0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image()

		z <- x0 %>% model@encoder()
		z <- z$mean()
		u <- z %>% model@decoder()

		mu_nfr <- u[, 1]
		sd_nfr <- u[, 2]
		sd_nfr <- sd_nfr + 1
		mu_mono <- u[, 3]
		sd_mono<- u[, 4]
		sd_mono <- sd_mono + 1

		browser()

		model
	}
)

