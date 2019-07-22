predict.vae <- function(model, x, batch_size = 2^14){

	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples
	sample_dim <- window_dim * num_samples

	batch_size2 <- round(batch_size / num_samples)	# per row
	starts <- seq(1, window_dim, by = batch_size2)
	ends <- starts + batch_size2 - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	if (model$prior == 'gmm'){

		Z <- NULL

		for (b in 1:n_batch){

			flog.info(sprintf('prediction | batch=%4.d/%4.d', b, n_batch))

			i <- starts[b]:ends[b]
			window_dim2 <- length(i)

			X <- mcols(x)$counts[i, , drop = FALSE] %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>% 
				tf$reshape(shape(window_dim2, num_samples, model$input_dim, model$feature_dim)) %>%
				tf$reshape(shape(window_dim2 * num_samples, model$input_dim, model$feature_dim)) %>%
				tf$expand_dims(axis = 3L)

			Y <- mcols(x)$coverage[i, , , drop = FALSE]	%>%
				tf$cast(tf$float32) %>%
				tf$transpose(c(0L, 2L, 1L)) %>%
				tf$reshape(shape(window_dim2 * num_samples, model$feature_dim)) %>%
				tf$expand_dims(axis = 2L)

			Z_x <- model$encoder$vplot(X)$loc
			Z_y <- model$encoder$coverage(Y)$loc
			Zb <- tf$concat(list(Z_x, Z_y), axis = 1L) %>%
				tf$reshape(shape(window_dim2, num_samples, 2 * model$latent_dim)) 
			Zb <- Zb %>% as.array()
			Zb <- aperm(Zb, c(1, 3, 2))
			Z <- abind(Z, Zb, along = 1)

			browser()
			i <- 1
			yy <- c(50, 200, 400, 600, 670)
			image(matrix(mcols(gr)$counts[i, , ], 32, 32), col = colorpanel(100, low = 'black', high = 'white'))
			x <- model$decoder$vplot(Zb[i, , drop = FALSE])$sample(1000L) %>% as.array()
			image(t(drop(colSums(x, dim = 1))), col = colorpanel(100, low = 'black', high = 'white'))
			plot(mcols(gr)$coverage[i, , 1], type = 'b')
			y <- model$decoder$coverage(Zb[i, , drop = FALSE])$sample(1000L) %>% as.array()
			plot(drop(colSums(y, dim = 1)), type = 'b')

			browser()
		}
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))

	Z
} # predict.vae
