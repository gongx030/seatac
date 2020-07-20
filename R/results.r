#' results
setMethod(
	'results',
	signature(
		x = 'Vplots',
		targets = 'GRangesList'
	),
  function(
		x,
		targets,
		...
	){

		browser()

		if (!all(granges(control) == granges(treatment)))
			stop('genomic ranges in control and treatment must be equal')

		nfr <- which(control@centers > 0 & control@centers <= 100)
		mono_nucleosome <-  which(control@centers > 180 & control@centers <= 247)

		n_control <- rowSums(control$counts)
		n_treatment <- rowSums(treatment$counts)

		nfr_control_mean <- rowMeans(control$distance_mean[, nfr])
		nfr_treatment_mean <- rowMeans(treatment$distance_mean[, nfr])

		nfr_mean <- nfr_treatment_mean - nfr_control_mean

		mono_control_mean <- rowMeans(control$distance_mean[, mono_nucleosome])
		mono_treatment_mean <- rowMeans(treatment$distance_mean[, mono_nucleosome])

		mono_mean <- mono_treatment_mean - mono_control_mean

		x <- granges(control)
		mcols(x) <- data.frame(
			window_id = 1:length(control),
			nfr_treatment = nfr_treatment_mean,
			nfr_control = nfr_control_mean,
			nfr_diff = nfr_mean,
			mono_treatment = mono_treatment_mean,
			mono_control = mono_control_mean,
			mono_diff = mono_mean
		)
		x
	}
)

