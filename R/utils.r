#' Specify a Guassian kernel
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
