setMethod(
	'cluster_vplot',
	signature(
		x = 'GRanges'
	),
	function(x, k = 3, nstart = 5, iter_max = 10, ...){

		nc <- ncol(x$counts)	# n_bins_per_window * n_intervals * n_samples

		# (n_bins_per_window * n_intervals * n_samples) ~ (n_intervals * n_intervals)
		C2G <- sparseMatrix(
			i = 1:nc,
			j = rep(1:(metadata(x)$n_bins_per_window * metadata(x)$n_intervals), times = metadata(x)$n_samples),
			dims = c(nc, metadata(x)$n_bins_per_window * metadata(x)$n_intervals)
		)

		# read counts: (n_bins_per_window * n_intervals) ~ n_samples
		W <- colSums(x$counts) %>% 
			matrix(nrow = metadata(x)$n_bins_per_window * metadata(x)$n_intervals, ncol = metadata(x)$n_samples)
		W <- W / rowSums(W)	# ratio of reads of each grid across all samples
		W[is.na(W)] <- 0

		r <- rep(1:(each = metadata(x)$n_bins_per_window * metadata(x)$n_intervals), metadata(x)$n_samples)

		# expected read count matrix: k ~ n_bins_per_window ~ n_intervals ~ n_samples
		E <- x$counts %*% C2G # length(x) ~ (n_intervals * n_bins_per_window)
		E <- E[, r] %*% Diagonal(x = c(W))

		X <- x$counts - E 	# raw deviations for each interval as the input
		XX <- rowSums(X^2)	

		res <- bplapply(seq_len(nstart), function(s){

			for (iter in seq_len(iter_max)){

				if (iter == 1){	# randomly initialize the centroids
					membership <- sample(1:k, length(x), replace = TRUE)
				}

				# a binary matrix for motif ~ cluster
				M <- sparseMatrix(
					i = 1:length(x),
					j = membership,
					dims = c(length(x), k) 
				) 
				M <- M %*% Diagonal(x = 1 / colSums(M))

				C <- t(M) %*% X	# k ~ (n_bins_per_window * n_intervals * n_samples)
				CC <- rowSums(C^2)

				D <- XX - 2 * X %*% t(C) + matrix(CC, nrow = nrow(X), ncol = k, byrow = TRUE)
				membership <- max.col(-D)

			}

			loss <- sum(D[cbind(1:nrow(X), membership)])	# the kmeans loss
			flog.info(sprintf('clustering | nstart=%2.d/%2.d | loss=%15.3f', s, nstart, loss))

			list(
				membership = membership, 
				loss = loss
			)
		})

		s <- which.min(sapply(res, function(y) y$loss))
		membership <- res[[s]]$membership

		x$deviations <- X

		x$membership <- membership

		metadata(x)$k <- k
		x
	}
)	


