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
load_pretrained_vplot_vae_model <- function(block_size = 240L, latent_dim = 10L, filters0 = 128L, min_reads = 15L, training = 'mix15'){

	n_intervals <- 48L

	message(sprintf('load_pretrained_vplot_vae_model | %s | block_size=%d', training, block_size))

	model_id <- sprintf('https://s3.msi.umn.edu/gongx030/projects/seatac/models/training=%s_block_size=%d_latent_dim=%d_filters0=%d_min_reads=%d', training, block_size, latent_dim, filters0, min_reads)
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
	res <- model(tf$random$uniform(shape(1L, n_intervals, model$n_bins_per_block, 1L)), tf$random$uniform(shape(1L, n_intervals)))
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


#'
#' @export
#'
VplotsList <- function(...){
	x <- list(...)
	new('VplotsList', x)
}

#'
#' @export
#'
setAs('ANY', 'VplotsList', function(from) {
	as(from, 'SimpleList')
})


#' Downsample V-plot
#'
#' Downsample a dense V-plot to a saprse V-plot, where the remaining number of reads is specified as `num_reads`.
#'
downsample_vplot <- function(x, num_reads = 1L){

	batch_size <- x$shape[[1]]

  j <- x %>%
		tf$reshape(shape(batch_size, -1L)) %>% 
		tf$math$log() %>%	 # to logit
		tf$random$categorical(num_reads) %>%
		tf$reshape(shape(batch_size * num_reads, 1L))

	i <- matrix(c(1:batch_size) - 1L, batch_size, num_reads) %>%
		tf$reshape(shape(batch_size * num_reads, 1L)) %>%
		tf$cast(tf$int64)

	y <- tf$SparseTensor(
		indices = tf$concat(list(i, j), axis = 1L),
		values = rep(1, batch_size * num_reads),
		dense_shape = c(batch_size, x$shape[[2]] * x$shape[[3]])
	) %>%
		tf$sparse$reorder() %>%
		tf$sparse$to_dense(validate_indices = FALSE) %>%	# not checking duplicated indices
		tf$reshape(shape(batch_size, x$shape[[2]], x$shape[[3]], 1L)) %>%
		scale_vplot()
	y

} # downsample_vplot


#' get_fragment_size
#'
get_fragment_size <- function(x){
	fragment_size <- x %>%
		tf$reduce_sum(shape(0L, 2L, 3L))

	fragment_size  <- (fragment_size / tf$reduce_sum(fragment_size)) %>%
		tf$expand_dims(0L) %>%
		tf$`repeat`(repeats = x$shape[[1]], axis = 0L)

	fragment_size	
} # get_fragment_size


