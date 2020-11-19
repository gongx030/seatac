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

		assay_names <- names(assays(x))
		assay_sum <- lapply(assay_names, function(i) tf$zeros(shape(length(classes), n_bins_per_block * x@n_intervals)))
		names(assay_sum) <- assay_names

		field_names <- colnames(rowData(x))
		if (length(field_names) > 0){
			field_sum <- lapply(field_names, function(h){
				dim_h <- dim(rowData(x)[[h]])
				if (dim_h[2] == n_blocks_per_window){
					if (length(dim_h) == 2)
						target_shape <- 1L
					else
						target_shape <- dim_h[-c(1:2)]
					tf$zeros(shape(length(classes), target_shape))
				}else if (dim_h[2] == x@window_size){
					tf$zeros(shape(length(classes), block_size))
				}else
					stop(sprintf('dim(rowData(x)$%s)[2] must be either %d or %d', h, x@window_size, n_blocks_per_window))
			})
			names(field_sum) <- field_names
		}

		freq <- tf$zeros(shape(length(classes)))	# number of motif hits

		for (i in 1:length(batches)){

			message(sprintf('summarize_vplot | batch=%6.d/%6.d', i, length(batches)))
			b <- batches[[i]]

			bins <- x[b] %>% 
				granges() %>%
				resize(fix = 'center', width = x@window_size - block_size + x@bin_size) %>%
				slidingWindows(x@bin_size, x@bin_size) %>%
				unlist()

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

			for (h in 1:length(assay_names)){

				BV <- assays(x[b])[[h]] %>%
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

				BV <- BV %>%
					tf$reshape(c(BV$shape[[1]], -1L))

				KV <- tf$sparse$sparse_dense_matmul(KB, BV)	# TF sites ~ Vplot
				CV <- tf$sparse$sparse_dense_matmul(CK, KV) # classes ~ Vplot

				assay_sum[[h]] <- assay_sum[[h]] + CV	# aggregated V-plot for each motif

			}

			if (length(field_names) > 0){

				for (h in 1:length(field_names)){
					dim_h <- dim(rowData(x)[[h]])
					if (dim_h[2] == n_blocks_per_window){

						BV <- rowData(x[b])[[h]] %>%
							as.array() %>%
							tf$cast(tf$float32) %>%
							tf$reshape(c(length(b) * n_blocks_per_window, -1L))
	
					}else if (dim_h[2] == x@window_size){
						BV <- rowData(x[b])[[h]] %>%
							as.matrix() %>%
							tf$cast(tf$float32)  %>%
							tf$expand_dims(-1L) %>%
							tf$expand_dims(-1L) %>%
							tf$image$extract_patches(
								sizes = c(1L, block_size, 1L, 1L),
								strides = c(1L, x@bin_size, 1L, 1L),
								rates = c(1L, 1L, 1L, 1L),
								padding = 'VALID'
							) %>%
							tf$squeeze(axis = 2L) %>%
							tf$reshape(c(length(b) * n_blocks_per_window, -1L))
					}

  	      KV <- tf$sparse$sparse_dense_matmul(KB, BV) # TF sites ~ Vplot
					CV <- tf$sparse$sparse_dense_matmul(CK, KV) # classes ~ Vplot
					field_sum[[h]] <- field_sum[[h]] + CV # aggregated V-plot for each motif
				}
			}
	
			freq <- freq + CK %>% tf$sparse$reduce_sum(1L) 	# number of motif hits

		}

		for (h in 1:length(assay_names)){
			assay_sum[[h]] <- assay_sum[[h]] / freq[h]
		}

		if (length(field_names) > 0){
			for (h in 1:length(field_names)){
				field_sum[[h]] <- field_sum[[h]] / freq[h]
			}
		}
		
		se <- SummarizedExperiment(
			assays = lapply(assay_sum, as.matrix)
		)

		if (length(field_names) > 0){
			for (h in field_names){
				SummarizedExperiment::rowData(se)[[h]] <- field_sum[[h]] %>% as.array()
			}
		}

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
