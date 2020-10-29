#' Summarize V-plot from a list of genomic range sets
#' 
#' Summarize V-plot from a list of genomic range sets
#'
#' @param x a Vplots object 
#' @param annotation a GRangesList object of motif binding sites
#' @param batch_size batch size
#' @param block_size block size
#'
#' @return a SummarizedVplots object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'summarize_vplot',
	signature(
		x = 'Vplots',
		annotation = 'GRangesList'
	), 
	function(
		x, 
		annotation,
		block_size,
		batch_size = 256L
	){

		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		classes <- names(annotation)
		annotation <- unlist(annotation)
		labels <- factor(names(annotation), classes)
		annotation <- annotation %>% 
			resize(width = 1L, fix = 'center')

		batches <- cut_data(length(x), batch_size)

		counts <- tf$zeros(shape(length(classes), n_bins_per_block * x@n_intervals))	# total observed counts of the V-plot for each TF
		freq <- tf$zeros(shape(length(classes)))	# number of motif hits

		for (i in 1:length(batches)){

			message(sprintf('summarize_vplot | batch=%6.d/%6.d', i, length(batches)))
			b <- batches[[i]]

			BV <- assays(x[b])$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					x@n_intervals,
					x@n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				extract_blocks_from_vplot(n_bins_per_block)	%>% # batch ~ n_blocks_per_window ~ n_intervals ~ n_bins_per_block ~ 1L
				tf$reshape(c(length(b) * n_blocks_per_window, x@n_intervals, n_bins_per_block, 1L))

			not_empty <- tf$reduce_sum(BV, shape(1L, 2L, 3L)) > 0	# has reads
			BV <- BV %>% tf$boolean_mask(not_empty)	# non-empty V-plot

			BV <- BV %>%
				tf$reshape(c(BV$shape[[1]], -1L))

			bins <- x[b] %>% 
				granges() %>%
				resize(fix = 'center', width = x@window_size - block_size + x@bin_size) %>%
				slidingWindows(x@bin_size, x@bin_size) %>%
				unlist()
			
			bins <- bins[not_empty %>% as.logical()]	# bins that have non-empty V-plot

			j <- annotation %over% bins
			n <- sum(j) # number of motif sites that have non-empty V-plot

			mm <- findOverlaps(annotation[j], bins) %>% 
				as.matrix()

			KB <- tf$SparseTensor(
				indices = mm - 1L,	# to zero-based
				values = rep(1, nrow(mm)),
				dense_shape = c(length(annotation[j]), length(bins))
			) # motif sites ~ bins

			CK <- tf$SparseTensor(
				indices = cbind(as.integer(labels)[j], 1:n) - 1L, # to zero-based
				values = rep(1, n),
				dense_shape = c(length(classes), n)
			) # classes ~ motif sites

			KV <- tf$sparse$sparse_dense_matmul(KB, BV)	# TF sites ~ Vplot
			CV <- tf$sparse$sparse_dense_matmul(CK, KV) # classes ~ Vplot

			counts <- counts + CV	# aggregated V-plot for each motif
			freq <- freq + CK %>% tf$sparse$reduce_sum(1L) 	# number of motif hits

		}

		se <- SummarizedExperiment(
			assays = list(
				counts = counts %>% as.matrix()
			)
		)
		SummarizedExperiment::rowData(se)$freq <- as.numeric(freq)
		SummarizedExperiment::rowData(se)$class <- classes

		new(
			'SummarizedVplots',
			se,
			fragment_size_range  = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			bin_size = x@bin_size,
			window_size = block_size,
			n_intervals = x@n_intervals,
			n_bins_per_window = n_bins_per_block,
			breaks = x@breaks,
			centers = x@centers,
			positions = seq(x@bin_size, block_size, by = x@bin_size) - (block_size / 2)
		)

	}
) # summarize_vplot
