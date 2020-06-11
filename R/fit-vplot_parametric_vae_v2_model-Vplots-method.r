#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_parametric_vae_v2_model',
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

		flog.info(sprintf('input windows:%d', length(x)))

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

				xc <- x[b]$counts %>%
					array_reshape(c(	
						length(b),
						metadata(x)$n_bins_per_window,
						metadata(x)$n_intervals
					))

				# read counts for each data point
				m <- xc@vals %>% 
					tf$cast(tf$float32) %>%
					tf$reshape(c(length(xc@vals), 1L))

				# a binary mapping from data point to V-plot (number of data points ~ batch_size)
				xc <- xc@subs
				g <- sparseMatrix(i = 1:nrow(xc), j = xc[, 1], dims = c(nrow(xc), length(b))) %>%
					as.matrix() %>%
					tf$cast(tf$float32)

				# read position relative to the center for each point
				y <- metadata(x)$centers[xc[, 3]] %>%
					tf$cast(tf$float32) %>%
					tf$reshape(c(nrow(xc), 1L))

				# binary index matrix of data point ~ window (relative to center) assignment
				h <- sparseMatrix(i = 1:nrow(xc), j = xc[, 2], dims = c(nrow(xc), metadata(x)$n_bins_per_window)) %>%
					as.matrix() %>%
					tf$cast(tf$float32)

				# weight for each v-plot
				w <- 1 / n[xc[, 1]]	%>%	
					tf$cast(tf$float32) %>%
					tf$reshape(c(nrow(xc), 1L))

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					posterior <- xi %>% 
						model@encoder()

					z <- posterior$sample()

					v <- z %>% 	# batch size ~ n_bins_per_window
						model@decoder()

					u <- tf$tensordot(g, v, 1L)

					y_pred <- (u * h) %>%
						tf$reduce_sum(1L, keepdims = TRUE)

					loss_likelihood <- (-w * m * (tfd_normal(y_pred, model@sigma0)$log_prob(y))) %>%
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

      if (epoch == 1 || epoch %% 10 == 0){
				plot(y[, 1], y_pred[, 1], main = epoch)
			}

			flog.info(sprintf('training %s | epoch=%4.d/%4.d | total_loss_likelihood=%13.7f | total_loss_kl=%13.7f | total_loss=%13.7f', class(model), epoch, epochs, total_loss_likelihood, total_loss_kl, total_loss))

		}

		model
	}
)

