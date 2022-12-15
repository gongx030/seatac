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


#' get_nucleosome_score
#' 
#' Calcualting the nucleosome score (eta) from Vplots
#'
#' @param x a tensorflow.tensor object 
#' @param beta0 Scale factor for calculating the nucleosome score 
#' @param beta1 Offset factor for calculating the nucleosome score 
#' 
#' @return a tensorflow.tensor object of the calculated nucleosome score 
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
get_nucleosome_score <- function(x, beta0 = -8.23990, beta1 = 13.19315){

	1 / (1 + 1 / exp(beta0 + beta1 * x))

} # get_nucleosome_score



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


#' block_center
#'
#' Label the center of a x-length sequence
#' @param x length of the sequence
#' @return a binary vector indicating the center of the sequence with length x
#' @export
#'
block_center <- function(x){
	y <- rep(FALSE, x)
	if (x %% 2 == 0){
		y[(x / 2):(x / 2 + 1)] <- TRUE
	}else
		y[(x + 1) / 2] <- TRUE
	y
}



#' get_fragment_size_per_batch
#'
#' Get the aggregated fragment size distribution per batch
#' @param x a Vplot object
#'
get_fragment_size_per_batch <- function(x){

	batch <- rowData(x)$batch %>%
		factor(x@samples) %>%
		as.numeric() %>%
		matrix(length(x), 1L) %>%
		tf$cast(tf$int32) %>%
		tf$math$subtract(1L) %>%
		tf$one_hot(x@n_samples) %>%
		tf$squeeze(1L)

	fragment_size <- assays(x)$counts %>%
		as.matrix() %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window)) %>%
		tf$reduce_sum(2L) # fragment size per sample

	fragment_size <- batch %>%
		tf$matmul(fragment_size, transpose_a = TRUE) %>%
		scale01()  # average fragment size per batch

	fragment_size
}


#' kernel_smoothing_1d
#' 
#' 1D kernel smoothing
#'
#' @param x a 2D tensor
#' @param size kernel size
#' @param sigma standard deviation of smoothing kernel
#'
#' @return a smoothed 2D tensor
#'
kernel_smoothing_1d <- function(x, size, sigma){

	d <- tfp$distributions$Normal(0, sigma)
	kernel <- d$prob(tf$range(start = -size, limit = size + 1L, dtype = tf$float32))
	kernel <- kernel / tf$reduce_sum(kernel)
	kernel <- kernel %>% tf$reshape(shape(-1L, 1L, 1L))
	x <- x %>%
		tf$expand_dims(2L) %>%
		tf$nn$conv1d(kernel, stride = 1L, padding = "SAME") %>%
		tf$squeeze(2L)
	x

}

#' standarize_1d
#' 
#' Standarize a 2D tensor
#'
#' @param x a 2D tensor
#'
#' @return a standarized 2D tensor
#'
standarize_1d <- function(x){
	x_mean <- x %>% tf$reduce_mean(1L, keepdims = TRUE)
	x_sd <- x %>% tf$math$reduce_std(1L, keepdims = TRUE)
	x <- (x - x_mean) / x_sd
	x
}

#' Execute code that traverses a dataset until an out of range condition occurs
#'
#' https://www.rdocumentation.org/packages/tfdatasets/versions/2.7.0/topics/until_out_of_range
#' As of 5/17/2022, `until_out_of_range` from the tfdatasets failed on colab, due to the following error:
#' Error: tensorflow.python.framework.errors_impl.OutOfRangeError: End of sequence [Op:IteratorGetNext]
#'
#' @param expr an expression
#'
until_out_of_range2 <- function(expr){
	break_error <- tryCatch(eval(parse(text = "break")), error = function(e) e)
  expr <- substitute(expr)
	envir <- parent.frame()
	tryCatch({
		while (TRUE) eval(expr, envir = envir)
	}, error = function(e) {
		NULL
#   if (!identical(e$message, break_error$message))
#     out_of_range_handler(e)
  })
}
