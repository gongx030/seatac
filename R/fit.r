
#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_predictor_model',
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

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		S <- x$sequence %>%
			as.character() %>%
			strsplit('')
		S <- do.call('rbind', S)
		S <- S %>% 
			factor(c('A', 'C', 'G', 'T')) %>%
			as.numeric() %>%
			matrix(nrow(S), ncol(S))
		S <- S - 1

		for (epoch in seq_len(epochs)){

			total_loss <- 0

			for (i in seq_len(n_batch)){

				b <- starts[i]:ends[i]

				xi <- S[b, ] %>%
					tf$cast(tf$int32)

				si <- (x$sample_id[b] - 1) %>%
					tf$cast(tf$int32)
					
				yi <- x[b]$counts  %>%
					as.matrix() %>%
					reticulate::array_reshape(c(		# convert into a C-style array
						length(b),
						metadata(x)$n_bins_per_window, 
						metadata(x)$n_intervals, 
						1L
					)) %>%
					tf$cast(tf$float32)

				with(tf$GradientTape(persistent = TRUE) %as% tape, {

					zi <- list(xi, si) %>% 
						model@encoder()

					yi_pred <- zi %>%
						model@decoder()

					loss <- (yi - yi_pred)^2 %>%
						tf$reduce_sum(axis = shape(1L, 2L, 3L)) %>%
						tf$reduce_mean()
				})

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

			flog.info(sprintf('training | epoch=%4.d/%4.d | total_loss=%9.1f', epoch, epochs, total_loss))

			if (epoch %% 5 == 0){

				ii <- sample(1:length(x), 1000)
				zii <- model %>% encode(x[ii], batch_size = batch_size)

				yii <- umap(zii)$layout
				plot(yii, pch = 21, bg = x$sample_id)

				yii_pred <- zii %>%
					tf$cast(tf$float32) %>%
					model@decoder()

				yii_pred %>%
					tf$reduce_sum(c(0L)) %>%
					tf$squeeze() %>%
					as.matrix() %>%
					t() %>%
					image(main = epoch)

				colSums(x[ii]$counts) %>%
					matrix(metadata(x)$n_bins_per_window, metadata(x)$n_intervals) %>%
					image()

			}
		}
	}
) # fit
