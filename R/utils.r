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
load_pretrained_vplot_vae_model <- function(block_size = 240L, latent_dim = 10L, filters0 = 128L, min_reads = 15L){

	n_intervals <- 48L

	message(sprintf('load_pretrained_vplot_vae_model | a pretrained model on GM12878 | block_size=%d', block_size))

	model_id <- sprintf('https://s3.msi.umn.edu/gongx030/projects/seatac/models/model=Vae_block_size=%d_latent_dim=%d_filters0=%d_min_reads=%d', block_size, latent_dim, filters0, min_reads)
	model_index_file <- sprintf('%s.index', model_id)
	model_data_file <- sprintf('%s.data-00000-of-00001', model_id)

	local_model_id <- tempfile()
	local_model_index_file <- sprintf('%s.index', local_model_id)
	local_model_data_file <- sprintf('%s.data-00000-of-00001', local_model_id)

	tryCatch({
		download.file(model_index_file, local_model_index_file)
	}, error = function(e){
		stop(sprintf('failed downloading %s', model_index_file))
	})

	tryCatch({
		download.file(model_data_file, local_model_data_file)
	}, error = function(e){
		stop(sprintf('failed downloading %s', model_data_file))
	})

	model <- VaeModel(block_size = block_size, n_intervals = n_intervals)
	res <- model(tf$random$uniform(shape(1L, n_intervals, model$n_bins_per_block, 1L)))
	load_model_weights_tf(model, local_model_id)
	new('VaeModel', model = model)

} # load_pretrained_vplot_vae_model


#' cut data
#' 
#' Cut a sequence into small batches
#' @param n Length of the sequence
#' @param batch_size Batch size
#' @return a list of sequence indices, where the length of the list is the number of batches.
#' @author Wuming Gong (gongx030@umn.edu)
#'
cut_data <- function(n, batch_size){
	starts <- seq(1, n, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > n] <- n
	n_batch <- length(starts)
	lapply(1:n_batch, function(i) starts[i]:ends[i])
}


#' compute_z_score
#'
compute_z_score <- function(x, pseudo_count = 1L){

	counts <- assays(x)$counts %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(length(x), x@n_intervals, x@n_bins_per_window))

	counts <- counts + pseudo_count

	# background probability along the fragment size dimension
	p <- tf$reduce_sum(counts, axis = shape(0L, 2L), keepdims = TRUE) / tf$reduce_sum(counts, keepdims = TRUE)
	n <- tf$reduce_sum(counts, 2L, keepdims = TRUE) # counts along the fragment size dimension per k-mer

	# Z-score per position per fragment size interval per k-mer
	Z <- (counts - n * p) / tf$math$sqrt(n * p * (1 - p)) # k-mer ~ fragment size ~ bins

	SummarizedExperiment::assays(x)$z <- Z %>% tf$reshape(shape(length(x), x@n_intervals * x@n_bins_per_window)) %>% as.matrix()

	x

} # compute_z_score

#' scale_vplot 
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
scale_vplot <- function(x){
	w <- x %>% tf$reduce_sum(shape(1L), keepdims = TRUE)  # sum of reads per bin
	x <- x / tf$where(w > 0, w, tf$ones_like(w))  # scale so that the sum of each bar is one (softmax)
	x
}


#' 
#' @export
#'
vplot2nucleosome <- function(x, is_nucleosome, is_nfr, scale = -10, offset = -0.95){
	di <- x %>%
		tf$boolean_mask(is_nucleosome, axis = 1L) %>%
		tf$reduce_sum(1L) %>%
		tf$squeeze(2L)

	nfr <- x %>%
		tf$boolean_mask(is_nfr, axis = 1L) %>%
		tf$reduce_sum(1L) %>%
		tf$squeeze(2L)

	1 / (1 + tf$math$exp(scale * (di / (di + nfr)  + offset)))

} # vplot2nucleosome

#'
#'
get_bsgenome <- function(x){
	if (x == 'hg19'){
		require(BSgenome.Hsapiens.UCSC.hg19)
		BSgenome.Hsapiens.UCSC.hg19
	}else if (x == 'mm10'){
		require(BSgenome.Mmusculus.UCSC.mm10)
		BSgenome.Mmusculus.UCSC.mm10
	}else
		spritnf('unknown genome:%s', x) %>% stop()
} # get_bsgenome

