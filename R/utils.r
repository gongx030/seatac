#' Specify a Guassian kernel
#' Adopted from https://stackoverflow.com/questions/52012657/how-to-make-a-2d-gaussian-filter-in-tensorflow
#'
gaussian_kernel <- function(size = 10, mu = 1, std = 1){
	# Make Gaussian Kernel with desired specs.

	d <- tfd_normal(mu, std)
	vals <- d$prob(tf$range(start = -size, limit = size + 1, dtype = tf$float32))
	gauss_kernel <- tf$einsum('i,j->ij', vals, vals)
	gauss_kernel <- gauss_kernel / tf$reduce_sum(gauss_kernel)

	# Expand dimensions of `gauss_kernel` for `tf.nn.conv2d` signature.
	gauss_kernel <- gauss_kernel[, , tf$newaxis, tf$newaxis]
	gauss_kernel
}


find_fragment_size_mixture <- function(x){

	x <- Reduce('c', x)

	fs <- colSums(x$counts) %>%
		matrix(x@n_bins_per_window, x@n_intervals) %>%
		colSums()

	d <- data.frame(
		fragment_size = rep(x@centers, fs)
	)

	m <- flexmix(
		fragment_size ~ 1,
		data = d,
		k = 2,
		model = list(FLXMRglm(family = 'Gamma'), FLXMRglm(family = 'gaussian'))
	)
	
	if (parameters(m, component = 1, model = 2)[1, ] > parameters(m, component = 2, model = 2)[1, ]){
		nfr_component <- 2
		mono_component <- 1
	}else{
		nfr_component <- 1
		mono_component <- 2
	}

	mu <- parameters(m, component = mono_component, model = 2)[1, ]
	sigma <- parameters(m, component = mono_component, model = 2)[2, ]
	shape <-  parameters(m, component = nfr_component, model = 1)[2, ]
	scale <-  1 / parameters(m, component = nfr_component, model = 1)[1, ] / shape

	list(mu = mu, sigma = sigma, shape = shape, scale = scale)
	
}

generate_random_vplots <- function(n, fsd, cd, n_bins_per_window){

	n_intervals <- length(fsd)

	# sampling total number of counts per vplot
	counts <- sample(1:length(cd), n, replace = TRUE, prob = cd)

	id <- rep(1:n, counts)

	# sampling fragment size for each read in each vplot
	fs <- sample(1:n_intervals, sum(counts), replace = TRUE, prob = fsd)

	# sampling window positions
	pos <- sample(1:n_bins_per_window, sum(counts), replace = TRUE)

	df <- data.frame(id = id, fragment_size = fs, position = pos) %>%
		group_by(id, fragment_size, position) %>%
		tally()

	x <- tf$SparseTensor(
		indices = cbind(df$id - 1, df$fragment_size - 1, df$position - 1), 
		values = df$n,
		dense_shape = c(as.integer(n), as.integer(n_intervals), as.integer(n_bins_per_window))
	) %>%
		tf$sparse$to_dense() %>%
		tf$cast(tf$float32) %>%
		tf$expand_dims(3L)

	x
}

extract_blocks_from_vplot <- function(x, n_bins_per_block){

	n_intervals <- x$shape[[2]]

	y <- x %>%
		tf$image$extract_patches(
			sizes = c(1L, n_intervals, n_bins_per_block, 1L),
			strides = c(1L, 1L, 1L, 1L),
			rates = c(1L, 1L, 1L, 1L),
			padding = 'VALID'
		) %>%
		tf$squeeze(axis = 1L)

	y <- y %>%
		tf$reshape(c(y$shape[[1]], y$shape[[2]], n_intervals, n_bins_per_block)) %>%
		tf$expand_dims(-1L)

	y
}

#' https://stackoverflow.com/questions/50706431/please-how-to-do-this-basic-thing-with-tensorflow
#'
reconstruct_vplot_from_blocks <- function(x){

	batch <- x$shape[[1]]
	n_blocks_per_window <- x$shape[[2]]
	n_intervals <- x$shape[[3]]
	n_bins_per_block <- x$shape[[4]]
	window_size <- as.integer(n_blocks_per_window + n_bins_per_block - 1)

	padding <- matrix(as.integer(c(0, 0, 0, 0, 0, n_blocks_per_window * n_intervals)), 3, 2, byrow = TRUE)
	zeros <- tf$zeros(c(batch, n_blocks_per_window, n_intervals))

	w <- rep(n_bins_per_block, window_size)
	w[1:n_bins_per_block] <- w[1:n_bins_per_block] - n_bins_per_block:1 + 1
	w[(window_size - n_bins_per_block + 1):window_size] <- w[(window_size - n_bins_per_block + 1):window_size] - 1:n_bins_per_block + 1

	w <- tf$constant(w, dtype = tf$float32) %>%
		tf$reshape(c(1L, 1L, window_size, 1L))

	y <- x %>% 
		tf$transpose(perm = c(0L, 1L, 3L, 2L, 4L)) %>%
		tf$reshape(c(batch, n_blocks_per_window, n_bins_per_block * n_intervals)) %>%
		tf$pad(padding, 'CONSTANT')

	y <- tf$concat(list(y, zeros), axis = 2L)
	y <- y %>% tf$reshape(c(batch, -1L))
	y <- y[, 1:(y$shape[[2]] - n_blocks_per_window * n_intervals)]
	y <- y %>% tf$reshape(c(batch, n_blocks_per_window, window_size + 1L, n_intervals))
	y <- y %>% tf$transpose(perm = c(0L, 1L, 3L, 2L))

	y <- y[, , , 1:window_size]

	y <- y %>% 
		tf$reduce_sum(axis = 1L) %>%
		tf$expand_dims(-1L)  %>%
		tf$multiply(1 / w)

	y

}

scale01 <- function(x){
	x <- (x - rowMins(x)) / (rowMaxs(x) - rowMins(x))
	x[is.na(x)] <- 0
	x
}

#' get_angles
#' Adopted from https://www.tensorflow.org/tutorials/text/transformer
#'
get_angles <- function(pos, i, d_model){
	angle_rates <- 1 / (10000^ ( (2 * i %/% 2) / d_model ))
  pos %*% angle_rates
}

positional_encoding <- function(position, d_model){
	angle_rads <- get_angles(
		matrix(0:(position - 1), position, 1),
		matrix(0:(d_model - 1), 1, d_model),
		d_model
	)
	even <- seq(1, ncol(angle_rads), by = 2)
	angle_rads[, even] <- sin(angle_rads[, even])
	odd <- seq(2, ncol(angle_rads), by = 2)
	angle_rads[, odd] <- cos(angle_rads[, odd])
	angle_rads

} # positional_encoding


create_look_ahead_mask <- function(size){
  mask <- 1 - tf$linalg$band_part(tf$ones(shape(size, size)), -1L, 0L)
  mask  # (seq_len, seq_len)
}

scale_vplot <- function(x){
	w <- tf$reduce_sum(x, c(1L, 3L), keepdims = TRUE)  # number reads for each bar
	w <- tf$where(w > 0, 1 / w, 1)  # compute the scale value
	x <- x * w  # scale so that the sum of reads per bar is one (or zero, if no reads is available in that bar)
	x
}


select_blocks <- function(x, block_size, min_reads, batch_size = 128, with_kmers = FALSE, with_nucleoatac = FALSE){

	kmers <- NULL
	nucleoatac <- NULL

	starts <- seq(1, length(x), by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > length(x)] <- length(x)
	n_batch <- length(starts)

	res <- lapply(seq_len(n_batch), function(j){
		h <- starts[j]:ends[j]
	  inputs <- x[h] %>% prepare_blocks(block_size, min_reads, with_kmers = with_kmers, with_nucleoatac = with_nucleoatac)
		flog.info(sprintf('select blocks | n_batch=%5.d/%5.d | selected=%5.d', j, n_batch, inputs$vplots$shape[[1]]))
		inputs
	})

	vplots <- tf$concat(lapply(res, function(xx) xx$vplots), axis = 0L)
	flog.info(sprintf('select blocks | total selected=%5.d', vplots$shape[[1]]))

	if (with_kmers){
		kmers <- tf$concat(lapply(res, function(xx) xx$kmers), axis = 0L)
	}

	if (with_nucleoatac){
		nucleoatac <- tf$concat(lapply(res, function(xx) xx$nucleoatac), axis = 0L)
	}

	list(vplots = vplots,  kmers = kmers)
} # 


#' prepare_blocks
#'
prepare_blocks <- function(x, block_size, min_reads = 0, with_kmers = FALSE, with_nucleoatac = FALSE){

	if (with_kmers && is.null(x$kmers))
		stop('x$kmers must be available if with_kmers=TRUE')

	if (with_nucleoatac && is.null(x$nucleoatac))
		stop('x$nucleoatac must be available if with_nucleoatac=TRUE')

	n_bins_per_block <- as.integer(block_size / x@bin_size)

	y <- x %>%
		prepare_vplot() %>%
		extract_blocks_from_vplot(n_bins_per_block) 
	
	y <- y %>%
		tf$reshape(c(y$shape[[1]] * y$shape[[2]], y$shape[[3]], y$shape[[4]], 1L))

	h <- y %>%
		tf$math$count_nonzero(c(1L, 3L), dtype = tf$float32) %>% # number of non-zero pixel per bar
		tf$math$count_nonzero(1L) %>% # number of bar that have non-zero pixels
		tf$cast(tf$float32)

	include <- h >= min_reads

	y <- y %>% 
		tf$boolean_mask(include, axis = 0L)

	if (with_nucleoatac){
		browser()
	}

	if (with_kmers){
		g <- x$kmers %>%
			tf$cast(tf$int32) %>%
			tf$expand_dims(-1L) %>%
			tf$expand_dims(-1L) %>%
			tf$image$extract_patches(
				sizes = c(1L, block_size, 1L, 1L),
				strides = c(1L, x@bin_size, 1L, 1L),
					rates = c(1L, 1L, 1L, 1L),
				padding = 'VALID'
			) %>%
			tf$squeeze(axis = 2L)

		g <- g %>%
			tf$reshape(c(g$shape[[1]] * g$shape[[2]], block_size))

		g <- g %>% tf$boolean_mask(include, axis = 0L)
	}else{
		g <- NULL
	}

	list(vplots = y, kmers = g)

} # prepare_blocks

setMethod(
	'prepare_vplot',
	signature(
		x = 'Vplots'
	),
	function(
		x
	){

		y <- x$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32)

		y

	}
)
