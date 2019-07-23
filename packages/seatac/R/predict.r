predict.vae <- function(model, x, window_size = 320, batch_size = 256, step_size = 50){

	expand <- metadata(x)$expand
	window_dim <- length(x)
	num_samples <- metadata(x)$num_samples

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	if (model$prior == 'gmm'){

		Z <- NULL

		for (b in 1:n_batch){

			flog.info(sprintf('prediction | batch=%4.d/%4.d', b, n_batch))

			i <- starts[b]:ends[b]
			xi <- makeData(x[i], window_size = window_size, train = FALSE, min_reads_per_window = 0)
			window_dim2 <- length(xi)

			X <- mcols(xi)$counts %>%
				as.matrix() %>%
				tf$cast(tf$float32) %>% 
				tf$reshape(shape(window_dim2, model$input_dim, model$feature_dim)) %>%
				tf$expand_dims(axis = 3L)

			Y <- mcols(xi)$coverage %>% 
				tf$cast(tf$float32) %>%
				tf$expand_dims(axis = 2L)

			Z_x <- model$encoder$vplot(X)$loc
			Z_y <- model$encoder$coverage(Y)$loc
			Z <- tf$concat(list(Z_x, Z_y), axis = 1L)

			browser()
			X <- model$decoder$vplot(Z)$sample(100L) %>% tf$reduce_sum(axis = 0L)
			Y <- model$decoder$coverage(Z)$sample(100L) %>% tf$reduce_sum(axis = 0L)


			i <- 1
			yy <- c(50, 200, 400, 600, 670)
			image(t(matrix(as.matrix(X[i, , , 1]), 32, 32)), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
			axis(2, at = (yy - 50) / (670 - 50), label = yy)
			axis(1, at = c(0, 0.5, 1))
			image(matrix(as.matrix(mcols(xi)$counts[i, , , 1]), 32, 32), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
			axis(2, at = (yy - 50) / (670 - 50), label = yy)
			axis(1, at = c(0, 0.5, 1))
			plot(Y[i, , 1], type = 'b', lwd = 2)
			plot(mcols(xi)$coverage[i, , 1], type = 'b', lwd = 2)

#			y <- model$decoder$coverage(Zb[i, , drop = FALSE])$sample(1000L) %>% as.array()
#			plot(drop(colSums(y, dim = 1)), type = 'b', lwd = 2)

			browser()
		}
	}else
		stop(sprintf('model$prior %s is not supported', model$prior))

	Z
} # predict.vae
