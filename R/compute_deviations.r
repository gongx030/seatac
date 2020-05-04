setMethod(
	'compute_deviations',
	signature(
		x = 'GRanges'
	),
	function(x, gc_group = 5, theta = 1, resampling = 100, burnin = 10, ...){

		flog.info(sprintf('gc group(gc_group):%d', gc_group))
		flog.info(sprintf('theta:%.3f', theta))
		flog.info(sprintf('resampling:%d', resampling))

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

		# a binary matrix for motif ~ membership
		M <- sparseMatrix(
			i = 1:length(x),
			j = x$membership,
			dims = c(length(x),  metadata(x)$k)
		) 

		r <- rep(1:(each = metadata(x)$n_bins_per_window * metadata(x)$n_intervals), metadata(x)$n_samples)

		# expected read count matrix: k ~ n_bins_per_window ~ n_intervals ~ n_samples
		ME <- t(M) %*% x$counts %*% C2G # k ~ (n_intervals * n_bins_per_window)
		ME <- (ME[, r] %*% Diagonal(x = c(W))) %>% as.matrix()
		dim(ME) <- c(metadata(x)$k, metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)

		# observed read count matrix: k~ n_bins_per_window ~ n_intervals ~ n_samples
		MX <- (t(M) %*% x$counts) %>% as.matrix()
		dim(MX) <- c(metadata(x)$k, metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)

		Y <- MX - ME 	# raw deviations
		Y <- smooth_vplot(Y, theta = theta)

		# split intervals into groups based on the GC content
		gc_cutoffs <- quantile(x$gc_content, seq(0, 1, length.out = gc_group + 1))
		groups <- as.numeric(cut(x$gc_content, gc_cutoffs, include.lowest = TRUE))

		B <- matrix(NA, length(x), resampling)
		for(i in seq_len(gc_group)){
			j <- which(groups == i)
			B[j, ] <- sample(j, resampling * length(j), replace = TRUE)
		}

		# mean deviation of resampled intervals
		Y_resample_mean <- Y_resample_var <- array(
			0, 
			dim = c(
				metadata(x)$k, 
				metadata(x)$n_bins_per_window, 
				metadata(x)$n_intervals, 
				metadata(x)$n_samples
			)
		)

		D_resample_mean  <- D_resample_var <- matrix(
			0, 
			metadata(x)$k, metadata(x)$n_samples,
			dimnames = list(NULL, metadata(x)$samples)
		)

		# compute the expected vplot for permutated motifs
		for (i in seq_len(resampling)){

			flog.info(sprintf('resampling | %d/%d', i, resampling))

			MBX <- (t(M[B[, i], ]) %*% x$counts) %>% as.matrix()
			MBE <- t(M[B[, i], ]) %*% x$counts %*% C2G
			MBE <- (MBE[, r] %*% Diagonal(x = c(W))) %>% as.matrix()	# 	k ~ (n_intervals * n_bins_per_window * n_samples)

			Y_resample <- (MBX - MBE) %>%
				array(dim = c(
					metadata(x)$k, 
					metadata(x)$n_bins_per_window,
					metadata(x)$n_intervals,
					metadata(x)$n_samples
				))
					
			# Welford's online algorithm for computing variance of Y_resample
			previous_mean <- Y_resample_mean
			Y_resample_mean <- Y_resample_mean + (Y_resample - Y_resample_mean) / i
			Y_resample_var <- Y_resample_var + (Y_resample - Y_resample_mean) * (Y_resample - previous_mean)

			# Z-score of resampled vplot
			Z_resample <- (Y_resample - Y_resample_mean) / sqrt(Y_resample_var / (i - 1))
			Z_resample[is.infinite(Z_resample) | is.na(Z_resample)] <- 0

			# distance of resampled v-plot
			D_resample <- Z_resample^2 %>%
				aperm(c(1, 4, 2, 3)) %>%
	  		rowSums(dim = 2) 

			# Welford's online algorithm for computing variance of D_resample
			previous_mean <- D_resample_mean
			D_resample_mean <- D_resample_mean + (D_resample - D_resample_mean) / i
			D_resample_var <- D_resample_var + (D_resample - D_resample_mean) * (D_resample - previous_mean)
		}

		# Z-score of input vplot
		Z <- (Y - Y_resample_mean) / sqrt(Y_resample_var / (resampling - 1))
		Z[is.infinite(Z) | is.na(Z)] <- 0

		D <- Z^2 %>%
			aperm(c(1, 4, 2, 3)) %>%
			rowSums(dim = 2) 

		P <- 1 - pnorm(D, mean = D_resample_mean, sd = sqrt(D_resample_var / (resampling - 1)))
		dimnames(P) <- list(NULL, metadata(x)$samples)
		P_adj <- apply(P, 2, p.adjust, method = 'BH')

		dim(Y) <- c(metadata(x)$k, metadata(x)$n_bins_per_window * metadata(x)$n_intervals * metadata(x)$n_samples)
		dim(Z) <- c(metadata(x)$k, metadata(x)$n_bins_per_window * metadata(x)$n_intervals * metadata(x)$n_samples)

		MX <- (t(M %*% Diagonal(x = 1 / colSums(M))) %*% x$counts) %>% as.matrix()

		se <- SummarizedExperiment(
			assays = list(deviations = Y, z = Z, counts = MX),
		)
		metadata(se) <- metadata(x)
		metadata(se)$resampling <- resampling
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

