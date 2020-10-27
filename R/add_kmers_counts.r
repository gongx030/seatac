# add_kmers_counts
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_kmers_counts',
	signature(
		x = 'VplotsKmers',
		y = 'SummarizedKmers'
	),
	function(x, y, threshold = 1e-3, pseudo_count = 1){

		Z <- assays(y)$z%>%	# z score
			tf$cast(tf$float32) %>%
			tf$reshape(shape(length(y), y@n_intervals, y@n_bins_per_window)) 

		# significant positions where the chance of occurence is higher than background
		A <- tf$where(Z > qnorm(1 - threshold), 1, 0)
		A <- A %>% tf$transpose(shape(0L, 2L, 1L)) # kmer ~ n_bins_per_window ~ n_intervals
		
		# mapping from a window to surrounding bins
		blocks <- granges(x) %>% 	# GRanges for every window bin
			slidingWindows(x@bin_size, x@bin_size) %>%
			unlist()

		gr <- granges(x) %>%	# GRanges for every k-mer
			slidingWindows(1L, 1L) %>% # to every k-mer
			unlist() 
		mcols(gr)$kmers <- c(t(rowData(x)$kmers))	# (zero-based) k-mer index of kmer

		gr_bins <- gr %>%	# GRanges for every bins surronding every k-mer (aka kmer bins)
			resize(fix = 'center', width = y@window_size) %>% # to k-mer centric V-plot
			slidingWindows(x@bin_size, x@bin_size) %>% # bin the position dimension
			unlist()

		mcols(gr_bins)$kmers <- rep(mcols(gr)$kmers, each = y@n_bins_per_window)
		mcols(gr_bins)$position <- rep(1:y@n_bins_per_window, length(gr)) - 1L	# zero-based position

		include <- (tf$reduce_sum(A, 2L) > 0) %>%	# if the kmer / bin combination has any significant reads
			tf$gather_nd(cbind(gr_bins$kmers, gr_bins$position)) %>%	# map to every kmer bin
			as.logical()
			
		gr_bins <- gr_bins[include] # include the kmer bins that are associated with significant reads

		# retrieve the fragment size distribution for each kmer bin
		S <- tf$gather_nd(A, cbind(gr_bins$kmers, gr_bins$position))	# length(gr_bins) ~ n_intervals

		mm <- findOverlaps(blocks, gr_bins) %>% 
			as.matrix() %>%
			tf$cast(tf$int64)
		mm <- mm - 1L	# to zero-based

		G <- tf$sparse$SparseTensor(mm, rep(1, nrow(mm)), dense_shape = shape(length(blocks), length(gr_bins)))

		a <- tf$sparse$sparse_dense_matmul(G, S) %>%
			tf$reshape(shape(length(x), x@n_bins_per_window, x@n_intervals)) %>%
			tf$reshape(shape(length(x), x@n_bins_per_window * x@n_intervals)) %>%
			as.matrix() %>% 
			as('dgCMatrix')

		SummarizedExperiment::assays(x)$kmers_counts <- a

		x

	}
)
