#' Warmup and cosine decay scheduler, the schedular used in the SimCLR paper
#' 
#' Adopted from https://github.com/google-research/simclr/blob/67b562524223e65c0414a15c3bd7ac007d28f629/tf2/model.py#L78
#'
#' @param base_learning_rate Baseline learning rate (default: 0.3)
#' @param num_examples Number of total training samples
#' @param batch_size Batch size (default : 256)
#' @param warmup_epochs Warmup epochs (default: 10)
#' @param epochs Total training epochs
#' @return a WarmUpAndCosineDecay object
#'
#' @export
#'
WarmUpAndCosineDecay <- function(
	base_learning_rate = 0.3,
	num_examples = NULL,
	batch_size = 256L,
	warmup_epochs = 10,
	epochs = 100L
){

	stopifnot(!is.null(num_examples))

	schedule <- reticulate::PyClass(
		classname = 'WarmUpAndCosineDecay',
		inherit = tf$keras$optimizers$schedules$LearningRateSchedule,
		defs = list(
			`__init__` = function(
				self, 
				base_learning_rate,
				num_examples,
				name = NULL
			){
				
				super()$`__init__`()

				# the default base_learning_rate by SinCLR is https://github.com/google-research/simclr/blob/6bf69ce127ae33e181e1a6c5777c84570cb5d147/tf2/run.py#L37
				self$base_learning_rate <- base_learning_rate
				self$num_examples <- num_examples

				return(NULL)

			},
			`__call__` = function(self, step){

				warmup_steps <- floor(self$num_examples * warmup_epochs / batch_size)
				total_steps <- floor(self$num_examples * epochs / batch_size) + 1L

				scaled_lr <- self$base_learning_rate * batch_size / 256
				learning_rate <- step / warmup_steps * scaled_lr 

				cosine_decay <- tf$keras$experimental$CosineDecay(scaled_lr, total_steps - warmup_steps)

				learning_rate <- tf$where(step < warmup_steps, learning_rate, cosine_decay(step - warmup_steps))
				learning_rate
															          
			}
		)
	)
	schedule(base_learning_rate, num_examples)
}
