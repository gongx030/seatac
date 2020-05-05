setMethod(
	'cluster_vplot',
	signature(
		x = 'GRanges'
	),
	function(x, n_neighbors = 50, ...){

		cls <- x$latent %>%
			nng(k = n_neighbors, mutual = TRUE) %>%
			cluster_louvain()	%>%
			membership()

		x$membership <- as.numeric(factor(cls))
		metadata(x)$num_clusters <- max(x$membership)

		x
	}
)	


