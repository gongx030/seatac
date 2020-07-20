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

	n_intervals <- x$shape[1]

	y <- x %>%
		tf$image$extract_patches(
			sizes = c(1L, n_intervals, n_bins_per_block, 1L),
			strides = c(1L, 1L, 1L, 1L),
			rates = c(1L, 1L, 1L, 1L),
			padding = 'VALID'
		) %>%
		tf$squeeze(axis = 1L)

	y <- y %>%
		tf$reshape(c(y$shape[0], y$shape[1], n_intervals, n_bins_per_block)) %>%
		tf$expand_dims(-1L)

	y
}

#' https://stackoverflow.com/questions/50706431/please-how-to-do-this-basic-thing-with-tensorflow
#'
reconstruct_vplot_from_blocks <- function(x){

	n_blocks_per_window <- x$shape[1]
	n_bins_per_block <- x$shape[3]
	window_size <- as.integer(n_blocks_per_window + n_bins_per_block - 1)
	n_intervals <- x$shape[2]
	batch <- x$shape[0]

	padding <- matrix(as.integer(c(0, 0, 0, 0, 0, n_blocks_per_window * n_intervals)), 3, 2, byrow = TRUE)
	zeros <- tf$zeros(c(batch, n_blocks_per_window, n_intervals))

	w <- c(1:n_bins_per_block, rep(n_bins_per_block, window_size - 2 * n_bins_per_block), n_bins_per_block:1)
	w <- tf$constant(w, dtype = tf$float32) %>%
		tf$reshape(c(1L, 1L, window_size, 1L))

	y <- x %>% 
		tf$transpose(perm = c(0L, 1L, 3L, 2L, 4L)) %>%
		tf$reshape(c(batch, n_blocks_per_window, n_bins_per_block * n_intervals)) %>%
		tf$pad(padding, 'CONSTANT')

	y <- tf$concat(list(y, zeros), axis = 2L)
	y <- y %>% tf$reshape(c(batch, -1L))
	y <- y[, 1:(y$shape[1] - n_blocks_per_window * n_intervals)]
	y <- y %>% tf$reshape(c(batch, n_blocks_per_window, window_size + 1L, n_intervals))
	y <- y %>% tf$transpose(perm = c(0L, 1L, 3L, 2L))

	y <- y[, , , 1:window_size]

	y <- y %>% 
		tf$reduce_sum(axis = 1L) %>%
		tf$expand_dims(-1L)  %>%
		tf$multiply(1 / w)

	y

}


