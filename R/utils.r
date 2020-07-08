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
