#' split_dataset
#' 
#' Split a tfdataset object into training and testing sets
#'
#' @param x a tfdataset object
#' @param test_size The ratio of the testing set (default 0.15)
#' @param batch_size Batch size (default: 64L)
#' @return a list that include a training and a testing dataset, where both of them are 
#'					tfdataset object
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


#' cut data
#' 
#' Cut a sequence into small batches
#'
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
#' Scale the Vplot so that the sum of reads at each position is one. 
#' 
#' @param x Vplots in tensorflow.tensor object ([batch, height, width, 1])
#' @return scaled Vplots in tensorflow.tensor object ([batch, height, width, 1])
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
scale_vplot <- function(x){
	w <- x %>% tf$reduce_sum(shape(1L), keepdims = TRUE)  # sum of reads per bin
	x <- x / tf$where(w > 0, w, tf$ones_like(w))  # scale so that the sum of each bar is one (softmax)
	x
}


#' vplot2nucleosome
#' 
#' Calcualting the nucleosome score (eta) from Vplots
#'
#' @param x Vplots in tensorflow.tensor object ([batch, height, width, 1])
#' @param is_nucleosome Logistic indicator for the nucleosome region in Vplots (e.g. between 180 and 247) ([height])
#' @param is_nfr Logistic indicator for the NFR region in Vplots (e.g. < 100) ([height])
#' @param scale Scale factor for calculating the nucleosome score (default: -10)
#' @param offset Offset factor for calculating the nucleosome score (default: -0.95)
#' 
#' @return a tensorflow.tensor object of the calculated nucleosome score ([batch, width])
#' @author Wuming Gong (gongx030@umn.edu)
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


#' VplotsList
#' 
#' Build a VplotsList from a list of Vplots
#'
#' @param ... Arguments passed to list()
#' @export
#'
VplotsList <- function(...){
	x <- list(...)
	new('VplotsList', x)
}

#' setAs
#'  
#' @export
#'
setAs('ANY', 'VplotsList', function(from) {
	as(from, 'SimpleList')
})


#' downsample_vplot
#'
#' Downsample a dense V-plot to a saprse V-plot, where the remaining number of reads is specified as `num_reads`.
#'
#' @param x a tensorflow.tensor object of Vplots
#' @param num_reads The target number of reads (default: 1L)
#' @return down-sampled Vplots in tensorflow.tensor object
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
