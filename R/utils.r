#' extract_blocks_from_vplot
#'
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

#' reconstruct_vplot_from_blocks
#'
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

#' scale01
#'
scale01 <- function(x){
	x_min <- tf$reduce_min(x, 1L, keepdims = TRUE)
	x_max <- tf$reduce_max(x, 1L, keepdims = TRUE)
	w <- x_max - x_min
	w <- tf$where(w > 0, w, tf$ones_like(w))	# scale so that the sum of each bar is one
	(x - x_min) / w
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


#' split_dataset
#'
split_dataset <- function(x, test_size = 0.15, batch_size = 64L){

	n <- as.numeric(x$cardinality()) # total number of samples
	n_train <- as.integer((1 - test_size) * n)  # train samples

	train <- x %>%
		dataset_take(n_train) %>%
		dataset_batch(batch_size)

	test <- x %>%
		dataset_skip(n_train) %>%
		dataset_batch(batch_size)

	list(train = train, test = test)

} # split_dataset


#' load_pretrained_vae_model
#'
#' Load pretrained sequence agnostic VAE model for V-plot
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
load_pretrained_vplot_vae_model <- function(n_intervals = 48L, block_size = 240L){
	flog.info('load_pretrained_vplot_vae_model | a pretrained model on GM12878')
	dir <- system.file('extdata', package = 'seatac')
	model <- VaeModel(block_size = block_size, n_intervals = n_intervals)
	res <- model(tf$random$uniform(shape(1L, n_intervals, model$n_bins_per_block, 1L)))
	load_model_weights_tf(model, sprintf('%s/model', dir))
	new('VaeModel', model = model)
} # load_pretrained_vplot_vae_model


#' 
#'
cut_data <- function(n, batch_size){
	starts <- seq(1, n, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > n] <- n
	n_batch <- length(starts)
	lapply(1:n_batch, function(i) starts[i]:ends[i])
}



