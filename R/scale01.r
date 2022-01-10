#' scale01
#' 
#' Scale the column a matrix to [0, 1]
#' @param x a dgCMatrix object
#' @return a scaled dgCMatrix object
#' @export
#'
setMethod(
	'scale01',
	signature(
		x = 'dgCMatrix'
	),
	function(
		x
	){
		x %>%
			as.matrix() %>%
			scale01() %>%
			as('dgCMatrix')

	}
)

#' scale01
#' 
#' Scale the column a matrix to [0, 1]
#' @param x matrix object
#' @return a scaled matrix object
#' @export
#'
setMethod(
	'scale01',
	signature(
		x = 'matrix'
	),
	function(
		x
	){
		rn <- rownames(x)
		cn <- colnames(x)
		x <- x %>% 
			tf$cast(tf$float32) %>%
			scale01() %>%
			as.matrix()
		rownames(x) <- rn
		colnames(x) <- cn
		x
	}
)

#' scale01
#' 
#' Scale the column a matrix to [0, 1]
#' @param x a tensorflow.tensor object
#' @return a scaled tensorflow.tensor object
#' @export
#'
setMethod(
	'scale01',
	signature(
		x = 'tensorflow.tensor'
	),
	function(
		x
	){
		x_min <- tf$reduce_min(x, 1L, keepdims = TRUE)
	  x_max <- tf$reduce_max(x, 1L, keepdims = TRUE)
		w <- x_max - x_min
		w <- tf$where(w > 0, w, tf$ones_like(w))  # scale so that the sum of each bar is one
		(x - x_min) / w
	}
)

#' scale01
#' 
#' Scale the column a matrix to [0, 1]
#' @param x a numeric vector
#' @return a scaled numeric vector
#' @export
#'
setMethod(
	'scale01',
	signature(
		x = 'numeric'
	),
	function(
		x
	){
		matrix(x, nrow = 1, ncol = length(x)) %>%
			scale01() %>%
			as.numeric()
	}
)

