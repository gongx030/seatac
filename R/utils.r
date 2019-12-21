to_sparse_tensor <- function(x){
	nr <- nrow(x)
	nc <- ncol(x)
	x <- as.matrix(summary(x))	# sparse matrix representation (i, j, x)
	tf$SparseTensor(indices = x[, 1:2] - 1, values = x[, 3], dense_shape = shape(nr, nc)) %>%
		tf$sparse$reorder() # https://www.tensorflow.org/api_docs/python/tf/sparse/SparseTensor
}
