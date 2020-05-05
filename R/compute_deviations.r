setMethod(
	'compute_deviations',
	signature(
		x = 'GRanges'
	),
	function(x, gc_group = 5, resampling = 100, ...){

		flog.info(sprintf('gc group(gc_group):%d', gc_group))
		flog.info(sprintf('resampling:%d', resampling))

		if (is.null(x$membership))
			stop('membership must be specified')

		if (is.null(x$latent))
			stop('latent must be specified')

		nc <- ncol(x$latent)	# latnet_dim * n_samples
		latent_dim <- nc / metadata(x)$n_samples

		# a binary matrix for length(x) ~ membership
		M <- sparseMatrix(
			i = 1:length(x),
			j = x$membership,
			dims = c(length(x),  metadata(x)$num_clusters)
		) 
		M <- M %*% Diagonal(x = 1 / colSums(M))

		# split intervals into groups based on the GC content
		gc_cutoffs <- quantile(x$gc_content, seq(0, 1, length.out = gc_group + 1))
		groups <- as.numeric(cut(x$gc_content, gc_cutoffs, include.lowest = TRUE))

		B <- matrix(NA, length(x), resampling)
		for(i in seq_len(gc_group)){
			j <- which(groups == i)
			B[j, ] <- sample(j, resampling * length(j), replace = TRUE)
		}

		Y <- t(M) %*% x$latent %>%
			as.matrix()

		D <- Y^2
		dim(D) <- c(metadata(x)$num_clusters, latent_dim, metadata(x)$n_samples)
		D <- D %>%
			aperm(c(1, 3, 2)) %>%
			rowSums(dim = 2)

		# mean deviation of resampled intervals
		Y_resample_mean <- Y_resample_var <- matrix(0, metadata(x)$num_clusters, nc)

		D_resample_mean  <- D_resample_var <- matrix(0, metadata(x)$num_clusters, metadata(x)$n_samples)

		# compute the expected vplot for permutated intervals
		for (i in seq_len(resampling)){

			flog.info(sprintf('resampling | %d/%d', i, resampling))

			Y_resample <- t(M[B[, i], ]) %*% x$latent %>%
				as.matrix()

			# Welford's online algorithm for computing variance of Y_resample
			previous_mean <- Y_resample_mean
			Y_resample_mean <- Y_resample_mean + (Y_resample - Y_resample_mean) / i
			Y_resample_var <- Y_resample_var + (Y_resample - Y_resample_mean) * (Y_resample - previous_mean)

			# distance of resampled v-plot
			D_resample <- Y_resample^2
			dim(D_resample) <- c(metadata(x)$num_clusters, latent_dim, metadata(x)$n_samples)
			D_resample <- D_resample %>%
				aperm(c(1, 3, 2)) %>%
				rowSums(dim = 2)

			# Welford's online algorithm for computing variance of D_resample
			previous_mean <- D_resample_mean
			D_resample_mean <- D_resample_mean + (D_resample - D_resample_mean) / i
			D_resample_var <- D_resample_var + (D_resample - D_resample_mean) * (D_resample - previous_mean)

		}

		P <- 1 - pnorm(D, mean = D_resample_mean, sd = sqrt(D_resample_var / (resampling - 1)))
		dimnames(P) <- list(NULL, metadata(x)$samples)
		P_adj <- apply(P, 2, p.adjust, method = 'BH')

		X_mean <- t(M) %*% x$counts 

		se <- SummarizedExperiment(
			assays = list(counts = t(M) %*% x$counts)
		)
		metadata(se) <- metadata(x)
		metadata(se)$resampling <- resampling
		rowData(se)$size <- table(x$membership)
		rowData(se)$p_value <- P
		rowData(se)$p_adj <- P_adj
		colData(se) <- expand.grid(
				position = metadata(x)$positions, 
				interval = metadata(x)$centers, 
				sample = metadata(x)$samples, 
				stringsAsFactors = FALSE
		) %>% DataFrame()
		se
	}
) # compute_deviations


#' smooth_vplot
#'
#' @param X input array
#' @return an array with the same dimensions as X
#' @author Wuming Gong (gongx030@umn.edu)
#'
smooth_vplot <- function(X, theta = 1){

	n_motifs <- dim(X)[1]
	n_bins_per_window <- dim(X)[2]
	n_intervals <- dim(X)[3]
	n_samples <- dim(X)[4]

	p <- expand.grid(
		 motif = seq_len(n_motifs),
		 sample = seq_len(n_samples)
	)
	Z <- bplapply(1:nrow(p), function(i){
		z <- X[p[i, 'motif'], , , p[i, 'sample']] %>%
			image.smooth(theta = theta) %>%
			pluck('z')
		z <- z / sum(z) * sum(X[p[i, 'motif'], , , p[i, 'sample']])
		z
	}) %>%
		unlist() %>%
		array(dim = c(
			n_bins_per_window,
			n_intervals,
			n_motifs,
			n_samples
    )) %>% 
		aperm(c(3, 1, 2, 4))
	Z
} # smooth_vplot

