#'
#' @export
#'
setMethod(
	'cluster_vplot',
	signature(
		x = 'Vplots'
	),
	function(
		x,
		k = 50,
		min_reads = 1
	){

		if (is.null(rowData(x)$latent))
			stop('rowData(x)$latent must be available')

		if (is.null(rowData(x)$n_reads))
			stop('rowData(x)$n_reads must be available')

		n_windows <- dim(rowData(x)$latent)[1]
		n_blocks_per_window <- dim(rowData(x)$latent)[2]

		z <- rowData(x)$latent %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(n_windows * n_blocks_per_window, dim(rowData(x)$latent)[3])) %>%
			as.matrix()

		n_reads <- rowData(x)$n_reads %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(n_windows * n_blocks_per_window)) %>%
			as.matrix()

		has_reads <- n_reads >= min_reads
		n <- sum(has_reads)

		message(sprintf('cluster_vplot | get.knn(k=%d)', k))
		gk <- get.knn(z[has_reads, ], k)

		A <- sparseMatrix(i = rep(1:n, k), j = c(gk$nn.index), x = 1 / (1 + c(gk$nn.dist)), dims = c(n, n))
		A <- (A + t(A)) / 2
		A <- graph.adjacency(A, mode = 'undirected', weighted = TRUE)

		message('cluster_vplot | cluster_louvain')
		g <- cluster_louvain(A)

		membership <- rep(-1, n_windows * n_blocks_per_window)
		membership[has_reads] <- g$membership

		SummarizedExperiment::rowData(x)$membership <- membership %>%
			tf$cast(tf$int32) %>%
			tf$reshape(shape(n_windows, n_blocks_per_window)) %>%
			as.matrix()
		x
	}
) # cluster_vplot

