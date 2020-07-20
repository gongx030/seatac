setMethod(
	'cluster',
	signature(
		x = 'VplotsFitted'
	),
	function(x, k, batch_size = 32L, ...){

		n <- length(x)

		z <- x$latent
		dim(z) <- c(n * x@model@n_blocks_per_window, x@model@latent_dim)

		km <- kmeans(z, k, ...)

		cls <- sparseMatrix(i = 1:nrow(z), j = km$cluster, dims = c(nrow(z), k)) %>%
			as('dgCMatrix') %>%
			as.matrix() %>%
			array(dim = c(n, x@model@n_blocks_per_window, k)) %>%
			tf$cast(tf$float32)

		starts <- seq(1, length(x), by = batch_size)
		ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		prototypes <- tf$zeros(c(k, x@n_intervals, x@model@n_bins_per_block, 1L))
		gr <- GRanges()

		for (i in 1:n_batch){

			b <- starts[i]:ends[i]

			xi <- x[b] %>%
				prepare_vplot() %>%
				extract_blocks_from_vplot(x@model@n_bins_per_block)  %>%
				tf$clip_by_value(clip_value_min = 0, clip_value_max = model@max_reads_per_pixel) %>%
				tf$reshape(c(length(b) * model@n_blocks_per_window, model@n_intervals, model@n_bins_per_block, 1L))

			ci <- cls[b, , ] %>%
				tf$reshape(c(length(b) * model@n_blocks_per_window, k))	
			
			prototypes <- prototypes + tf$tensordot(tf$transpose(ci), xi, axes = 1L)

			gri <- granges(x[b]) %>% 
				slidingWindows(width = x@model@block_size, step = x@bin_size) %>%
				unlist()
			gr <- c(gr, gri)

			browser()

			if (i == 1)
				mem <- ci
			else
				mem <- list(mem, ci) %>% tf$concat(0L)

		}

		gr$cluster <- tf$argmax(mem, 1L) %>% as.matrix() %>% c()
		gr$cluster <- gr$cluster + 1


		w <- tf$linalg$diag(1 / tf$reduce_sum(cls, c(0L, 1L)))
		prototypes <- tf$tensordot(w, prototypes, axes = 1L)

		prototypes <- prototypes %>%
			tf$squeeze(-1L) %>%
			as.array()

		prototypes <- aperm(prototypes, c(1, 3, 2))

		new(
			'VplotsCluster',
			gr,
			fragment_size_range  = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			bin_size = x@bin_size,
			window_size = x@model@block_size,
			n_intervals = x@n_intervals,
			n_bins_per_window = x@model@n_bins_per_block,
			breaks = x@breaks,
			centers = x@centers,
			positions = seq(x@bin_size, x@model@block_size, by = x@bin_size) - (x@model@block_size / 2),
			prototypes = prototypes,
			k = k 
		)

	}
)

