#' load_model
#'
#' Load a pretrianed model
#' @param model an empty VaeModel
#' @param filename  Filename of the pretrained model
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'load_model',
	signature(
		model = 'VaeModel',
		filename = 'character'
	),
	function(
		model,
		filename
	){

		model_index_file <- sprintf('%s.index', filename)
		model_data_file <- sprintf('%s.data-00000-of-00001', filename)

		stopifnot(file.exists(model_index_file))
		stopifnot(file.exists(model_data_file))

		res <- model@model(
			tf$random$uniform(shape(1L, model@model$n_intervals, model@model$n_bins_per_block, 1L)), 
			tf$random$uniform(shape(1L, model@model$n_intervals))
		)

		load_model_weights_tf(model@model, filename)
		model
	}
)
