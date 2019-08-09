#' segment
#' Separate the vplot based on whether the center is NFR or nucleosome

segment <- function(x, method = 'z', cutoff = 0.05, batch_size = 1000){

	if (is.null(mcols(x)$fitted_counts))
		stop('fitted_counts must be available')

	k <- 3	# two states: nfr or nucleosome
	num_steps <- metadata(x)$n_bins_per_window
  feature_dim <- ncol(mcols(x)$fitted_counts) / num_steps
	window_dim <- length(x)
	input_dim <- metadata(x)$n_bins_per_window
	feature_dim <- metadata(x)$n_intervals
	bin_size <- metadata(x)$bin_size

	breaks <- seq(metadata(x)$fragment_size_range[1], metadata(x)$fragment_size_range[2], by = metadata(x)$fragment_size_interval)

	is_nfr <- breaks[-1] > 0 & breaks[-1] <= 100
	is_nucleosome <- breaks[-1] > 180 
	is_mono_nucleosome <- breaks[-1] >= 180 & breaks[-1] <= 247
	is_di_nucleosome <- breaks[-1] >= 315 & breaks[-1] <= 473
	is_tri_nucleosome <- breaks[-1] >= 558 & breaks[-1] <= 615

	X <- mcols(x)$fitted_counts
	dim(X) <- c(window_dim, input_dim, feature_dim)

	R <- log(rowMeans(X[, , is_mono_nucleosome | is_di_nucleosome | is_tri_nucleosome], dims = 2) + 0.1) - log(rowMeans(X[, , is_nfr], dims = 2) + 0.1)

	if (method == 'hmm'){

		starts <- seq(1, window_dim, by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > window_dim] <- window_dim
		n_batch <- length(starts)

		n_batch <- round(window_dim / batch_size)
		groups <- sample(1:n_batch, window_dim, replace = TRUE)
		S <- matrix(NA, nrow = nrow(X), ncol = ncol(X))

		for (b in 1:n_batch){
			flog.info(sprintf('segmenting | batch=%d/%d', b, n_batch))
			i <- groups == b
			ntimes <- rep(input_dim, sum(i))
			y <- c(t(R[i, , drop = FALSE]))
			m <- depmix(y ~ 1, nstates = k, family = gaussian(), ntimes = ntimes)
			fm <- depmixS4::fit(m, em = em.control(maxit = 100))
			state <- posterior(fm)[, 1]
			state <- order(sapply(1:k, function(j) mean(y[state == j])))[state]
			S[i, ] <- matrix(state, sum(i), input_dim, byrow = TRUE)
		}
		mcols(x)$state <- S

	}else if (method == 'z'){

		mcols(x)$nucleosome <- rowMeans(X[, , is_mono_nucleosome | is_di_nucleosome | is_tri_nucleosome], dims = 2)

		mcols(x)$z_score <- t((t(R) - colMeans(R)) / colSds(R))
		P <- pnorm(mcols(x)$z_score)

		S <- matrix(0, nrow = nrow(X), ncol = ncol(X))	# unknown
		S[mcols(x)$z_score  > 1] <-  1	# nucleosome
		S[mcols(x)$z_score < -1] <- -1	 # NFR
		mcols(x)$state <- S
	}
	x

} # segment


call_nucleosome <- function(x, center = 20){

	y <- Reduce('c', lapply(1:length(x), function(i){
		y <- rle(mcols(x)$state[i, ])
		starts <- cumsum(c(1, y$lengths))
		starts <- starts[-length(starts)]
		ends <- starts + y$lengths - 1
		n <- length(y$values)
		xi <- rep(x[i], n)
		xi <- narrow(xi, start = starts, end = ends)
		xii <- granges(xi)
		mcols(xii)$group <- mcols(xi)$group
		mcols(xii)$segment_id_per_window <- 1:n
		mcols(xii)$num_segments_per_window <- n
		mcols(xii)$state <- y$values
		mcols(xii)$window_id <- i
		xii
	}))
	mcols(y)$state <- factor(mcols(y)$state, c('nfr', 'nucleosome'))
	y

} # call_nucleosome
