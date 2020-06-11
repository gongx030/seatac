#' results
setMethod(
	'results',
	signature(
		control = 'Vplots',
		treatment = 'Vplots'
	),
  function(
		control,
		treatment,
		...
	){

		if (!all(granges(control) == granges(treatment)))
			stop('genomic ranges in control and treatment must be equal')

		nfr <- which(control@centers > 0 & control@centers <= 100)
		mono_nucleosome <-  which(control@centers > 180 & control@centers <= 247)

		n_control <- rowSums(control$counts)
		n_treatment <- rowSums(treatment$counts)

		nfr_control_mean <- rowMeans(control$distance_mean[, nfr])
		nfr_treatment_mean <- rowMeans(treatment$distance_mean[, nfr])

		nfr_mean <- nfr_treatment_mean - nfr_control_mean
		nfr_sd <- rowMeans(sqrt(control$distance_sd[, nfr]^2 + treatment$distance_sd[, nfr]^2))

		nfr_p <- pnorm(nfr_mean, mean = 0, sd = nfr_sd)
		nfr_p <- rowMins(cbind(nfr_p, 1 - nfr_p))

		mono_control_mean <- rowMeans(control$distance_mean[, mono_nucleosome])
		mono_treatment_mean <- rowMeans(treatment$distance_mean[, mono_nucleosome])

		mono_mean <- mono_treatment_mean - mono_control_mean
		mono_sd <- rowMeans(sqrt(control$distance_sd[, mono_nucleosome]^2 + treatment$distance_sd[, mono_nucleosome]^2))

		mono_p <- pnorm(mono_mean, mean = 0, sd = mono_sd)
		mono_p <- rowMins(cbind(mono_p, 1 - mono_p))

		x <- granges(control)
		mcols(x) <- data.frame(
			window_id = 1:length(control),
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

