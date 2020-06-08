
setMethod(
	'results',
	signature(
		x = 'GRanges',
		treatment = 'integer',
		control = 'integer'
	),
  function(
		x,
		treatment,
		control,
		...
	){

		if (treatment > metadata(x)$n_samples || control > metadata(x)$n_samples)
			stop(sprintf('treatment or control must be no greater than %d (number of samples)', metadata(x)$n_samples))

		nfr <- which(metadata(x)$centers > 0 & metadata(x)$centers <= 100)
		mono_nucleosome <-  which(metadata(x)$centers > 180 & metadata(x)$centers <= 247)

		i <- x$sample_id == control
		j <- x$sample_id == treatment

		n_control <- rowSums(x[i]$counts)
		n_treatment <- rowSums(x[j]$counts)

		nfr_control_mean <- rowMeans(x[i]$distance_mean[, nfr])
		nfr_treatment_mean <- rowMeans(x[j]$distance_mean[, nfr])

		nfr_mean <- nfr_treatment_mean - nfr_control_mean
		nfr_sd <- rowMeans(sqrt(x[i]$distance_sd[, nfr]^2 + x[j]$distance_sd[, nfr]^2))

		nfr_p <- pnorm(nfr_mean, mean = 0, sd = nfr_sd)
		nfr_p <- rowMins(cbind(nfr_p, 1 - nfr_p))

		mono_control_mean <- rowMeans(x[i]$distance_mean[, mono_nucleosome])
		mono_treatment_mean <- rowMeans(x[j]$distance_mean[, mono_nucleosome])

		mono_mean <- mono_treatment_mean - mono_control_mean
		mono_sd <- rowMeans(sqrt(x[i]$distance_sd[, mono_nucleosome]^2 + x[j]$distance_sd[, mono_nucleosome]^2))

		mono_p <- pnorm(mono_mean, mean = 0, sd = mono_sd)
		mono_p <- rowMins(cbind(mono_p, 1 - mono_p))

		x <- x[x$sample_id == 1]
		x <- granges(x)

		mcols(x) <- data.frame(
			window_id = 1:sum(i),
			nfr_treatment = nfr_treatment_mean,
			nfr_control = nfr_control_mean,
			nfr_diff = nfr_mean,
			nfr_pvalue = nfr_p,
			nfr_fdr = p.adjust(nfr_p, 'fdr'),
			mono_treatment = mono_treatment_mean,
			mono_control = mono_control_mean,
			mono_diff = mono_mean,
			mono_pvalue = mono_p,
			mono_fdr = p.adjust(mono_p, 'fdr')
		)
		x
	}
)

