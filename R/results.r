#' results
#'
#' Test the difference between two Vplots
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param type Tesing type ('phase', 'nucleosome', or 'vplots')
#' @param contrast this argument species a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)
#' @param ... Other arguments
#'
#' @return a GRanges object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'results',
	signature(
		model = 'VaeModel',
		x = 'Vplots'
	),
	function(
		model, 
		x,
		type = NULL,
		contrast = NULL,
		...
	){

		validate(model, x)

		stopifnot(!is.null(rowData(x)[['vae_z_mean']]))
		stopifnot(!is.null(rowData(x)[['vae_z_stddev']]))

		stopifnot(!is.null(type))
		stopifnot(is.character(contrast))
		stopifnot(length(contrast) == 3)

		field <- contrast[1]
		control <- contrast[2]
		treatment <- contrast[3]

		stopifnot(!is.null(x@dimdata[['sample']][[field]]))
		stopifnot(control %in% x@dimdata[['sample']][[field]])
		stopifnot(treatment %in% x@dimdata[['sample']][[field]])

		if (type == 'phase'){
			res <- results_phase(model, x, contrast, ...)
		}else if (type == 'nucleosome'){
			res <- results_nucleosome(model, x, contrast, ...)
		}else if (type == 'vplots'){
			stopifnot(!is.null(rowData(x)[['predicted_nucleosome']]))
			res <- results_vplots(x, field, treatment, control, ...)
		}else
			stop(sprintf('unknown type: %s', type))

		res

	}
)


#' results_phase
#'
#' Test the difference of the phase between two Vplots
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param contrast this argument species a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)
#' @param width Width of the centeral regions to be considered (default: 100L)
#' @param repeats Number of repeated sampling (default: 100L)
#' @param batch_size Batch size (default: 16L)
#'
#' @return a GRanges object
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
results_phase <- function(model, x, contrast, width = 100L, repeats = 100L, batch_size = 16L){

	field <- contrast[1]
	control <- contrast[2]
	treatment <- contrast[3]

	x <- x[rowData(x)[[field]] %in% c(treatment, control)]

  d <- model %>%
	  prepare_data(x) %>%
		tensor_slices_dataset() %>%
		dataset_batch(batch_size)
	iter <- d %>% make_iterator_one_shot()

	is_center <- x@positions >= -width / 2 & x@positions <= width /2
	fs_mean <- NULL
	fs_stddev <- NULL

	i <- 1
	res <- until_out_of_range({
		batch <- iterator_get_next(iter)
		batch$vplots <- batch$vplots %>%
			tf$sparse$to_dense() %>%
			scale01()
		res <- model@model(batch, training = FALSE)
		z <- res$posterior$sample(repeats) %>%
			tf$reshape(shape(repeats * res$z$shape[[1]], -1L)) 
		b <- tf$zeros(shape(z$shape[[1]]), dtype = tf$int64) %>% 
			tf$one_hot(model@model$n_batches)
		fs <- list(z, b) %>%
			tf$concat(1L) %>%
			model@model$decoder(training = FALSE) %>%
			tf$keras$activations$softmax(1L) %>%
			tf$boolean_mask(is_center, 2L) %>%
			tf$reduce_sum(2L) %>%
			tf$squeeze(2L) %>%
			tf$reshape(shape(repeats, res$z$shape[[1]], -1L))

		fs_mean <- c(fs_mean, fs %>% tf$reduce_mean(0L))
		fs_stddev <- c(fs_stddev, fs %>% tf$math$reduce_std(0L))
		if (i %% 10 == 0)
			sprintf('results_phase | batch=%6.d/%6.d', i, ceiling(nrow(x) / batch_size)) %>% message()
		i <- i + 1
	})

	fs_mean <- fs_mean %>% tf$concat(axis = 0L)
	fs_stddev <- fs_stddev %>% tf$concat(axis = 0L)

	rowData(x)[['fragment_size_mean']] <- as.matrix(fs_mean)
	rowData(x)[['fragment_size_stddev']] <- as.matrix(fs_stddev)

	k <- ncol(rowData(x)[['fragment_size_mean']])
	fs_control <- rowData(x[rowData(x)[[field]] == control])[['fragment_size_mean']] 
	fs_control_stddev <- rowData(x[rowData(x)[[field]] == control])[['fragment_size_stddev']] 
	fs_treatment <- rowData(x[rowData(x)[[field]] == treatment])[['fragment_size_mean']] 
	fs_treatment_stddev <- rowData(x[rowData(x)[[field]] == treatment])[['fragment_size_stddev']] 

	h <- ((fs_treatment - fs_control)^2 * (1 / (fs_control_stddev^2 + fs_treatment_stddev^2 ))) %>%
		rowSums()

	res <- granges(x[rowData(x)[[field]] == control])
	mcols(res) <- mcols(res)[c('id')]
	res$stat <- h
	res$pvalue <- 1 - pchisq(h, df = k)
	res$padj <- p.adjust(res$pvalue, method = 'BH')
	res

}

#' results_vplots
#'
#' Test the difference of Vplots by a Chi-squared test
#'
#' @param x a Vplots object
#' @param field The field in x@dimdata[['sample']] that defines the sample groups
#' @param treatment the name of the numerator level for the fold change
#' @param control the name of the denominator level for the fold change
#' @param width Width of the centeral regions to be considered (default: 100L)
#' @importFrom abind abind
#' @importFrom stats p.adjust
#'
#' @return a GRanges object
#'
results_vplots <- function(x, field, treatment, control, width = 100L){

	stopifnot(is.numeric(width) && width >= 0 && width <= x@window_size)

	latent_dim <- dim(rowData(x)[['vae_z_mean']])[2]

	i <- x@dimdata[['sample']][[field]] == treatment
	z_treatment <- rowData(x)[['vae_z_mean']][, i, , drop = FALSE] 
	z_treatment <- aperm(z_treatment, c(2L, 1L, 3L))
	z_treatment <- colMeans(z_treatment)
	z_stddev_treatment <- rowData(x)[['vae_z_stddev']][, i, , drop = FALSE] 
  z_stddev_treatment  <- aperm(z_stddev_treatment, c(2L, 1L, 3L))
  z_stddev_treatment <- colSums(z_stddev_treatment^2) %>% sqrt()
	nucleosome_treatment <- rowData(x)[['predicted_nucleosome']][, i, , drop = FALSE]
	nucleosome_treatment <- aperm(nucleosome_treatment, c(2L, 1L, 3L))
	nucleosome_treatment <- colMeans(nucleosome_treatment)

	i <- x@dimdata[['sample']][[field]] == control
	z_control <- rowData(x)[['vae_z_mean']][, i, , drop = FALSE] 
	z_control <- aperm(z_control, c(2L, 1L, 3L))
	z_control <- colMeans(z_control)
	z_stddev_control <- rowData(x)[['vae_z_stddev']][, i, , drop = FALSE] 
  z_stddev_control <- aperm(z_stddev_control, c(2L, 1L, 3L))
  z_stddev_control <- colSums(z_stddev_control^2) %>% sqrt()
	nucleosome_control <- rowData(x)[['predicted_nucleosome']][, i, , drop = FALSE]
	nucleosome_control <- aperm(nucleosome_control, c(2L, 1L, 3L))
	nucleosome_control <- colMeans(nucleosome_control)

	h <- ((z_treatment - z_control)^2 / (z_stddev_treatment^2 + z_stddev_control^2)) %>% 
		rowSums()

	pvalue_z <- 1 - pchisq(h, df = latent_dim)
	res <- granges(x)
	mcols(res) <- NULL
	res$pvalue_z <- pvalue_z
	res$padj <- p.adjust(pvalue_z)

	is_center <- x@dimdata$bin$position >= -width / 2 & x@dimdata$bin$position <= width / 2
	res$nucleosome_treatment <- rowMeans(nucleosome_treatment[, is_center, drop = FALSE])
	res$nucleosome_control <- rowMeans(nucleosome_control[, is_center, drop = FALSE])
	res$log_ratio <- log(res$nucleosome_treatment + .Machine$double.eps) - log(res$nucleosome_control + .Machine$double.eps)

	res

}

#' results_nucleosome
#'
#' Testing the difference of nucleosomes between two samples
#'
#' @param model a trained VaeModel object
#' @param x a Vplots object
#' @param contrast this argument species a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)
#' @param fragment_size_threshold Fragment size threshold for nucleosome reads
#' @param batch_size batch size (default: 16L)
#'
#' @return a GRanges object
#'
#' @export
results_nucleosome <- function(model, x, contrast, fragment_size_threshold = 150L, batch_size = 16L){

	field <- contrast[1]
	control <- contrast[2]
	treatment <- contrast[3]

	x <- x[rowData(x)[[field]] %in% c(treatment, control)]

  d <- model %>%
	  prepare_data(x) %>%
		tensor_slices_dataset() %>%
		dataset_batch(batch_size)
	iter <- d %>% make_iterator_one_shot()

	is_nucleosome <- x@centers >= fragment_size_threshold
	z_mean <- NULL
	z_stddev <- NULL
	nuc_nfr <- NULL

	i <- 1
	res <- until_out_of_range({

		batch <- iterator_get_next(iter)
		batch$vplots <- batch$vplots %>%
			tf$sparse$to_dense() %>%
			scale01()
		res <- model@model(batch, training = FALSE)

		z <- res$posterior$mean()

		z_mean <- c(z_mean, z)
		z_stddev <- c(z_stddev, res$posterior$stddev())

		b <- tf$zeros(shape(z$shape[[1]]), dtype = tf$int64) %>% 
			tf$one_hot(model@model$n_batches)

		y <- list(z, b) %>%
			tf$concat(1L) %>%
			model@model$decoder(training = FALSE) %>%
			tf$keras$activations$softmax(1L)

		nuc <- y %>%
			tf$boolean_mask(is_nucleosome, 1L) %>%
			tf$reduce_sum(1L)

		nfr <- y %>%
			tf$boolean_mask(!is_nucleosome, 1L) %>%
			tf$reduce_sum(1L)

		nn <- log((nuc + 1e-3) / (nfr + 1e-3)) %>%
			tf$squeeze(2L)

		nuc_nfr <- c(nuc_nfr, nn)

		if (i %% 10 == 0)
			sprintf('results_nucleosome | batch=%6.d/%6.d', i, ceiling(nrow(x) / batch_size)) %>% message()
		i <- i + 1
	})


	z_mean <- z_mean %>% tf$concat(axis = 0L)
	z_stddev <- z_stddev %>% tf$concat(axis = 0L)
	nuc_nfr <- nuc_nfr %>% tf$concat(axis = 0L)

	rowData(x)[['z_mean']] <- as.matrix(z_mean)
	rowData(x)[['z_stddev']] <- as.matrix(z_stddev)
	rowData(x)[['nuc_nfr']] <- as.matrix(nuc_nfr)

	nuc_nfr_control <- rowData(x[rowData(x)[[field]] == control])[['nuc_nfr']] 
	nuc_nfr_treatment <- rowData(x[rowData(x)[[field]] == treatment])[['nuc_nfr']] 

	z_control <- rowData(x[rowData(x)[[field]] == control])[['z_mean']] 
	z_control_stddev <- rowData(x[rowData(x)[[field]] == control])[['z_stddev']] 
	z_treatment <- rowData(x[rowData(x)[[field]] == treatment])[['z_mean']] 
	z_treatment_stddev <- rowData(x[rowData(x)[[field]] == treatment])[['z_stddev']] 

	h <- ((z_treatment - z_control)^2 / (z_control_stddev^2 + z_treatment_stddev^2)) %>%
		rowSums()
	pvalue_z <- 1 - pchisq(h, df = model@model$latent_dim)

	res <- granges(x[rowData(x)[[field]] == control])
	mcols(res) <- mcols(res)[c('id')]
	res$pvalue_z <- pvalue_z

	res$nucleosome_diff <- 2 / (1 + exp(-1 * (nuc_nfr_treatment - nuc_nfr_control))) - 1
	res

}
