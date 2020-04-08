setGeneric('compute_deviations', function(x, ...) standardGeneric('compute_deviations'))

setMethod(
	'compute_deviations',
	signature(
		x = 'GRanges'
	),
	function(x, gc_group = 5, resampling = 100, ...){

		nc <- ncol(x$counts)	# n_bins_per_window * n_intervals * n_samples

		# (n_bins_per_window * n_intervals * n_samples) ~ (n_intervals * n_intervals)
		C2G <- sparseMatrix(
			i = 1:nc,
			j = rep(1:(metadata(x)$n_bins_per_window * metadata(x)$n_intervals), times = metadata(x)$n_samples),
			dims = c(nc, metadata(x)$n_bins_per_window * metadata(x)$n_intervals)
		)

		# total reads per motif/grid combination
		X <- x$counts	%*% C2G # read counts: motifs ~ grid ( n_intervals * n_bins_per_window )

		# read counts: (n_bins_per_window * n_intervals) ~ n_samples
		W <- colSums(x$counts) %>% 
			matrix(nrow = metadata(x)$n_bins_per_window * metadata(x)$n_intervals, ncol = metadata(x)$n_samples)

		W <- W / rowSums(W)	# ratio of reads of each grid across all samples
		W[is.na(W)] <- 0

		# expected read counts (length(x) ~ n_intervals * n_bins_per_window * n_samples)
		E <- do.call('cbind', lapply(seq_len(metadata(x)$n_samples), function(i) X %*% Diagonal(x = W[, i])))

		# a binary matrix for motif ~ TF assignment
		M <- sparseMatrix(
			i = 1:length(x),
			j = as.numeric(x$tf_name),
			dims = c(length(x),  metadata(x)$n_motifs)
		) 

		# expected read count matrix: n_motifs ~ n_bins_per_window ~ n_intervals ~ n_samples
		ME <- (t(M) %*% E) %>% as.matrix()
		dim(ME) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)

		ME <- smooth_vplot(ME)

		# observed read count matrix: n_motifs ~ n_bins_per_window ~ n_intervals ~ n_samples
		MX <- (t(M) %*% x$counts) %>% as.matrix()
		dim(MX) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)

		Y <- MX - ME 	# raw deviations

		# split motif intervals into groups based on the GC content
		gc_cutoffs <- quantile(x$gc_content, seq(0, 1, length.out = gc_group + 1))
		groups <- as.numeric(cut(x$gc_content, gc_cutoffs, include.lowest = TRUE))

		# a permutation matrix (n_motifs ~ resampling)
		B <- matrix(NA, length(x), resampling)
		for (i in seq_len(gc_group)){
			j <- which(groups == i)
			B[j, ] <- sample(j, resampling * sum(groups == i), replace = TRUE)
		}

		# compute the expected vplot for permutated motifs
		MBE <- Reduce('+', bplapply(seq_len(resampling), function(i) (t(M[B[, i], ]) %*% E) %>% as.matrix()))  / resampling
		dim(MBE) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)
		MBE <- smooth_vplot(MBE)

		Y_resample <- bplapply(seq_len(resampling), function(i){
			MBX <- (t(M[B[, i], ]) %*% x$counts) %>% as.matrix()
			dim(MBX) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)
			MBX - MBE
		}) %>%	
			unlist()%>%
			array(dim = c(
				metadata(x)$n_motifs, 
				metadata(x)$n_bins_per_window,
				metadata(x)$n_intervals,
				metadata(x)$n_samples, 
				resampling
			))

		Y_resample_mean <- rowMeans(Y_resample, dims = 4)
		Y_resample_mean <- array(
			c(Y_resample_mean), 
			dim = c(
				metadata(x)$n_motifs, 
				metadata(x)$n_bins_per_window,
				metadata(x)$n_intervals,
				metadata(x)$n_samples, 
				resampling
			)
		)

		Y_resample_sd <- sqrt(rowSums((Y_resample - Y_resample_mean)^2, dims = 4) / resampling)

		Z <- (Y - rowMeans(Y_resample_mean, dims = 4)) / Y_resample_sd 

		dim(Z) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window * metadata(x)$n_intervals * metadata(x)$n_samples)
		dim(Y) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window * metadata(x)$n_intervals * metadata(x)$n_samples)

		positions <- seq(metadata(x)$bin_size, metadata(x)$window_size, by = metadata(x)$bin_size) - (metadata(x)$window_size / 2)

		se <- SummarizedExperiment(
			assays = list(deviations = Y, z = Z),
			rowData = data.frame(motif = metadata(x)$motifs)
		)
		metadata(se) <- metadata(x)
		rownames(se) <- metadata(x)$motifs
		metadata(se)$resampling <- resampling
		colData(se) <- expand.grid(
				position = positions, 
				interval = metadata(x)$centers, 
				sample = metadata(x)$samples, 
				stringsAsFactors = FALSE
		) %>% DataFrame()
		se
	}
)

#' smooth_vplot
#'
#' @param X
#' @return an array with the same dimensions as X
#' @author Wuming Gong (gongx030@umn.edu)
#'
smooth_vplot <- function(X){

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
			image.smooth(theta = 2) %>%
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

