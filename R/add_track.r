#' add_track
#'
#' @param object a bigwig file
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_track',
	signature(
		x = 'Vplots',
		object = 'character'
	),
	function(x, object, label){

		if (!file.exists(object))
			stop(sprintf('%s does not exist', object))

		cvg <- tryCatch({
			rtracklayer::import(object, which = reduce(x@rowRanges), as = 'RleList')
		}, error = function(e){
			stop(sprintf('%s cannot be imported by rtracklayer::import', objecct))
		})

		y <- cvg[x@rowRanges] %>% as.matrix()
		SummarizedExperiment::rowData(x)[[label]] <- y
		x
	}
) # add_track


#' add_track
#'
#' @param object a GRange object
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_track',
	signature(
		x = 'Vplots',
		object = 'GRanges'
	),
	function(x, object, label){

		cvg <- coverage(object)
		y <- cvg[x@rowRanges] %>% as.matrix()
		y[y > 0] <- 1
		y <- as(y, 'dgCMatrix')
		SummarizedExperiment::rowData(x)[[label]] <- y
		x
	}
) # add_track

