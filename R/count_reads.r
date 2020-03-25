setGeneric('count_reads', function(x, ...) standardGeneric('count_reads'))

#' count_reads
#'
setMethod(
	'count_reads',
	signature(
		x = 'GRanges'
	),
	function(x){

		nfr <- 1:100
		mono <- 180:247
		di <- 316:473
		tri <- 558:615

		x$counts <- cbind(
			nfr = rowSums(x$fragment_size_profile[, nfr]),
			mono_nucleosome = rowSums(x$fragment_size_profile[, mono]),
			di_nucleosome = rowSums(x$fragment_size_profile[, di]),
			tri_nucleosome = rowSums(x$fragment_size_profile[, tri])
		)
		x
	}
) # count_reads

