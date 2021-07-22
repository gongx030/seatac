#' prepare_data
#'
#' Prepare dataset for training and a V-plot model
#' 
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a Vplots object
#' @param weight Whether or not include positional weight
#'
#' @return a list that include `vplots`, `weight` and `batch`
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'prepare_data',
	signature(
		model = 'TransformerEncoderModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		...
	){

		d <- list()

		vplots <- assays(x)$counts %>%
			as.matrix() %>%
			tf$cast(tf$float32) %>%
			tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
			scale_vplot()

		coverage <- rowData(x)$coverage %>%
			tf$cast(tf$float32) 

		list(vplots = vplots, coverage = coverage)
	}
)


#' fit
#'
#' Fit a TransformerEncoderModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param warmup Warmup epochs (default: 50L)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#'
#' @export
#' @return a VaeModel
#'
setMethod(
	'fit',
	signature(
		model = 'TransformerEncoderModel',
		x = 'tf_dataset'
	),
	function(
		 model,
		 x,
		 batch_size = 128L,
		 epochs = 100L,
		 learning_rate = 1e-3,
		 compile = TRUE
	 ){

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		vplots_reconstrution_loss <- tf$keras$losses$BinaryCrossentropy(reduction = 'none')	# loss for the V-plot

		train_step <- function(x){
			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				res <- model@model(x, training = TRUE)
				loss_reconstruction_vplots <- vplots_reconstrution_loss(x$vplots, res$vplots) %>%
					tf$reduce_sum(shape(1L, 2L)) %>%
					tf$reduce_mean()
				loss <- loss_reconstruction_vplots
	 		})
			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss = loss,
				loss_reconstruction_vplots = loss_reconstruction_vplots
			)
		}


		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		for (epoch in seq_len(epochs)){
			loss <- NULL
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch)
				loss <- rbind(loss, sapply(res, as.numeric))
			})

			if (epoch == 1 || epoch %% 10 == 0){
				loss <- colMeans(loss)
				sprintf('epoch=%6.d/%6.d | %s', epoch, epochs, paste(sapply(1:length(loss), function(i) sprintf('%s=%13.7f', names(loss)[i], loss[i])), collapse = ' | ')) %>%
					message()
			}
		}
		model
	}
)

#' predict
#'
#' Predict counts
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param batch_size Batch size (default: 256L)
#'
#' @return a Vplots object with the mean and stddev of the latent representation
#'  (reducedDim(x, 'vae_z_mean') and reducedDim(x, 'vae_z_stddev'))
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'predict',
	signature(
		model = 'TransformerEncoderModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 256L, # v-plot per batch
		nfr_threshold = 100,
		...
	){

		d <- model %>%
			prepare_data(x, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()
		is_nfr <- x@centers <= nfr_threshold

		predicted_vplots <- NULL
		nucleosome_signal <- NULL

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch, training = FALSE)
			predicted_vplots <- c(predicted_vplots, res$vplots)
			nuc <- 1 - res$vplots %>% tf$boolean_mask(is_nfr, 1L) %>% tf$reduce_sum(1L) %>% tf$squeeze(2L)
#			s <- 1 / (1 + tf$math$exp(scale * (nfr) + offset))
#			s <- kernel_smoothing_1d(s, size, sigma)
#			s <- scale01(s)
			nucleosome_signal <- c(nucleosome_signal, nuc)
		})

		predicted_vplots <- predicted_vplots %>% 
			tf$concat(axis = 0L) %>%
			tf$reshape(shape(nrow(x), -1L)) %>% 
			as.matrix()
		dimnames(predicted_vplots) <- dimnames(x)
		assays(x)$predicted_counts <- predicted_vplots

		nucleosome_signal <- nucleosome_signal %>% tf$concat(axis = 0L)
		rowData(x)[['nucleosome_signal']] <- as.matrix(nucleosome_signal)

		x
	}
)

