#' build_classifier
#'
build_classifier <- function(x, core_width = 100){

	core_bin <- round(metadata(x)$n_bins_per_window / 2 - core_width / metadata(x)$bin_size):round(metadata(x)$n_bins_per_window / 2 + core_width / metadata(x)$bin_size)
	mono_region <- matrix(FALSE, metadata(x)$n_bins_per_window, metadata(x)$n_intervals)
	mono_region[core_bin, which(metadata(x)$mono_nucleosome)] <- TRUE
	nfr_region <- matrix(FALSE, metadata(x)$n_bins_per_window, metadata(x)$n_intervals)
	nfr_region[core_bin, which(metadata(x)$nfr)] <- TRUE
	nfr <- rowMeans(mcols(x)$counts[, which(nfr_region)])
	nfr <- nfr / mean(nfr)
	mono <- rowMeans(mcols(x)$counts[, which(mono_region)])
	mono <- mono / mean(mono)
	log_mono2nfr <- log(mono + 1) - log(nfr + 1)
	pos <- log_mono2nfr > quantile(log_mono2nfr, 0.75, na.rm = TRUE)
	neg <- log_mono2nfr < quantile(log_mono2nfr, 0.25, na.rm = TRUE)
	d <- data.frame(y = rep(c(1, 0), c(sum(pos), sum(neg))), rbind(mcols(x)$latent[pos, ], mcols(x)$latent[neg, ]))
	classifier <- glm(y ~ ., d, family = 'binomial')
	classifier	

} # build_classifier
