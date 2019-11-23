#' sample_windows 
#'
sample_windows <- function(x, window_size = 640, min_reads_per_window = 20, epochs = 50, batch_size = 128, steps_per_epoch = 20){

	flog.info(sprintf('training window size(window_size):%d', window_size))
	flog.info(sprintf('minimum PE reads per window(min_reads_per_window):%d', min_reads_per_window))
	flog.info(sprintf('training epochs(epochs):%d', epochs))
	flog.info(sprintf('batch size(batch_size):%d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch):%d', steps_per_epoch))

	bin_size <- metadata(x)$bin_size
	n_bins_per_window <- window_size / bin_size

	bin_starts <- 1:(metadata(x)$n_bins_per_window - n_bins_per_window + 1)
	seq_starts <- seq(1, metadata(x)$window_size - window_size + 1, by = metadata(x)$bin_size)
	n <- length(bin_starts)

	N <- epochs * batch_size * steps_per_epoch	# total number of training windows

	flog.info(sprintf('total number of training windows:%d', N))

	gs <- table(factor(sample(1:n, N, replace = TRUE), 1:n))	# number of training windows per start

	y <- lapply(1:n, function(i){

		flog.info('sampling training windows %3.d/%3.d(n=%6.d)', i, n, gs[i])
		pos <- rep(FALSE, metadata(x)$n_bins_per_window)
		pos[bin_starts[i]:(bin_starts[i] + n_bins_per_window - 1)] <- TRUE
		pos <- rep(pos, metadata(x)$n_intervals)
		counts <- mcols(x)$counts[, pos]
		num_reads <- Matrix::rowSums(counts)

		j <- sample(which(num_reads > min_reads_per_window), gs[i], replace = TRUE)

		list(
			counts = counts[j, , drop = FALSE], 
			seqnames = seqnames(x[j]),
			start = start(x[j]) + seq_starts[i] - 1,
			window_id = mcols(x[j])$window_id,
			group = mcols(x[j])$group
		)
	})

	gr <- GRanges(
		seqnames = Reduce('c', lapply(y, function(z) z$seqnames)), 
		range = IRanges(start = unlist(lapply(y, function(z) z$start)), width = window_size)
	)
	mcols(gr)$counts <- Reduce('rbind', lapply(y, function(z) z$counts))
	mcols(gr)$num_reads <- Matrix::rowSums(mcols(gr)$counts)
	mcols(gr)$window_id <- unlist(lapply(y, function(z) z$window_id))
	mcols(gr)$group <- unlist(lapply(y, function(z) z$group))

	gr <- gr[sample(1:length(gr))]	# shuffling

	# assigning batch index and step index
	mcols(gr)$epoch <- rep(1:epochs, each = batch_size * steps_per_epoch)
	mcols(gr)$step <- rep(1:steps_per_epoch, each = batch_size, times = epochs)

	metadata(gr) <- metadata(x)

	metadata(gr)$training <- list(
		epochs = epochs, 
		batch_size = batch_size, 
		steps_per_epoch = steps_per_epoch, 
		min_reads_per_window = min_reads_per_window
	)
	metadata(gr)$window_size <- window_size
	metadata(gr)$n_bins_per_window <- n_bins_per_window

	gr

} # sample_windows
