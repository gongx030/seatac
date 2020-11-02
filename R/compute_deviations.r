#' Compute deviations
#' 
#' Compute the V-plot deviations
#'
#' @param x a Vplots object 
#' @param annotation a GRangesList object of motif binding sites
#' @param model a pretrained VaeModel
#' @param batch_size_window batch size for processing windows
#' @param batch_size_block batch size for running VAE prediction
#' @param background number of background V-plots
#' @param permutation number of permutations
#'
#' @return a SummarizedVplots object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'compute_deviations',
	signature(
		x = 'Vplots',
		annotation = 'GRangesList',
		model = 'VaeModel'
	), 
	function(
		x, 
		annotation,
		model,
		batch_size_window = 128L,
		batch_size_block = 256L,
		background = 1000L,
		permutation = 100L
	){

		block_size <- model@model$block_size 

		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		classes <- names(annotation)
		annotation <- unlist(annotation)
		labels <- factor(names(annotation), classes)
		annotation <- annotation %>% 
			resize(width = 1L, fix = 'center')

		batches <- cut_data(length(x), batch_size_window)

		counts <- tf$zeros(shape(length(classes), n_bins_per_block * x@n_intervals))	# total observed counts of the V-plot for each TF
		freq <- tf$zeros(shape(length(classes)))	# number of motif hits

		Z <- tf$zeros(shape(length(classes), model@model$encoder$latent_dim))	# latent representation for each motif
		Zp <- tf$zeros(shape(permutation, length(classes), model@model$encoder$latent_dim))	# permutated latent representation 

		# pre-compute how many background V-plot should be sampled per batch
		bg <- sample(1:length(batches), background, replace = TRUE) %>%
			factor(1:length(batches)) %>%
			table()

		Z_bg <- list()

		for (i in 1:length(batches)){

			message(sprintf('compute_deviations | batch=%6.d/%6.d', i, length(batches)))

			b <- batches[[i]]

			BV <- assays(x[b])$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					x@n_intervals,
					x@n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				extract_blocks_from_vplot(n_bins_per_block)	%>% # batch ~ n_blocks_per_window ~ n_intervals ~ n_bins_per_block ~ 1L
				tf$reshape(c(length(b) * n_blocks_per_window, x@n_intervals, n_bins_per_block, 1L))

			not_empty <- tf$reduce_sum(BV, shape(1L, 2L, 3L)) > 0	# has reads

			BV <- BV %>% tf$boolean_mask(not_empty)	# non-empty V-plot

			BZ <- model %>% encode(BV, batch_size = batch_size_block)	# map to the latent space

			bins <- x[b] %>% 
				granges() %>%
				resize(fix = 'center', width = x@window_size - block_size + x@bin_size) %>%
				slidingWindows(x@bin_size, x@bin_size) %>%
				unlist()
			
			bins <- bins[not_empty %>% as.logical()]	# bins that have non-empty V-plot

			j <- annotation %over% bins
			n <- sum(j) # number of motif sites that have non-empty V-plot

			mm <- findOverlaps(annotation[j], bins) %>% 
				as.matrix()

			KB <- tf$SparseTensor(
				indices = mm - 1L,	# to zero-based
				values = rep(1, nrow(mm)),
				dense_shape = c(length(annotation[j]), length(bins))
			) # motif sites ~ bins

			CK <- tf$SparseTensor(
				indices = cbind(as.integer(labels)[j], 1:n) - 1L, # to zero-based
				values = rep(1, n),
				dense_shape = c(length(classes), n)
			) # classes ~ motif sites

			KZ <- tf$sparse$sparse_dense_matmul(KB, BZ)	# TF sites ~ latent
			CZ <- tf$sparse$sparse_dense_matmul(CK, KZ) # classes ~ latent

			Z <- Z + CZ	# aggregated latent for each motif

			BV <- BV %>%
				tf$reshape(c(BV$shape[[1]], -1L))

			KV <- tf$sparse$sparse_dense_matmul(KB, BV)	# TF sites ~ Vplot
			CV <- tf$sparse$sparse_dense_matmul(CK, KV) # classes ~ Vplot

			counts <- counts + CV	# aggregated V-plot for each motif

			# randomly sample background latent representation
			if (bg[i] > 0){
				m <- seq_len(BV$shape[[1]])%in% sample(seq_len(BV$shape[[1]]), bg[i], replace = TRUE)
				m <- tf$cast(m, tf$bool)
				Z_bg[[i]] <- BZ %>% tf$boolean_mask(m)
			}

			freq <- freq + CK %>% tf$sparse$reduce_sum(1L) 	# number of motif hits

			# permutation
			CK_perm <- tf$SparseTensor(
				indices = cbind(
					rep(1:permutation, each = n),
					rep(as.integer(labels)[j], permutation) %>% sample(),	# randomly permute the motif ~ bin relatinoship
					rep(1:n, permutation)
				) - 1L,
			 values = tf$ones(permutation * n),
			 dense_shape = shape(permutation, length(classes), n)
		 ) %>% # classes ~ permutation ~ TF sites
			tf$sparse$reshape(shape(permutation *  length(classes), n))

			CZ_perm <- tf$sparse$sparse_dense_matmul(CK_perm, KZ) %>%
				tf$reshape(shape(permutation, length(classes), -1L))

			Zp <- Zp + CZ_perm
		}

		Z_bg <- Z_bg[!sapply(Z_bg, is.null)]
		Z_bg <- Z_bg %>% tf$concat(0L)	# background latent representation

		Z_mean <- Z_bg %>% tf$reduce_mean(0L, keepdims = TRUE)	# mean latent vector
		Z_std <- Z_bg %>% tf$math$reduce_std(0L, keepdims = TRUE)	# std latent vector

		p <- tfp$distributions$MultivariateNormalDiag(loc = Z_mean, scale_diag = Z_std)

		Z <- Z / (freq %>% tf$expand_dims(1L))	# average latent for each motif
		Zp <- Zp / (freq %>% tf$expand_dims(0L) %>% tf$expand_dims(2L))	# average permuted latent for each motif

		# probability of mean latent vector under a multivariate normal distribution 
		d <- p$prob(Z) %>% tf$expand_dims(0L)

		# probability of permuted latent vector under a multivariate normal distribution 
		d_null <- p$prob(Zp)

		n <- (d - d_null > 0) %>% 
			tf$cast(tf$float32) %>% 
			tf$reduce_sum(0L)

		pvalue <- (n + 1) / permutation
		pvalue <- tf$math$minimum(pvalue, 1)

		X <- counts %>%
			tf$reshape(shape(length(classes), x@n_intervals, n_bins_per_block, 1L)) 
		res <- model %>% predict(X)

		predicted_counts <- res$vplots %>%
			tf$reshape(shape(length(classes), x@n_intervals * n_bins_per_block))

		se <- SummarizedExperiment(
			assays = list(
				counts = counts %>% as.matrix(),
				predicted_counts = predicted_counts %>% as.matrix()
			)
		)

		SummarizedExperiment::rowData(se)$freq <- as.numeric(freq)
		SummarizedExperiment::rowData(se)$class <- classes
		SummarizedExperiment::rowData(se)$pvalue <- as.numeric(pvalue)
		SummarizedExperiment::rowData(se)$p_adj <-  p.adjust(as.numeric(pvalue), method = 'BH')
		SummarizedExperiment::rowData(se)$nucleosome <- as.matrix(res$nucleosome)

		new(
			'SummarizedVplots',
			se,
			fragment_size_range  = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			bin_size = x@bin_size,
			window_size = block_size,
			n_intervals = x@n_intervals,
			n_bins_per_window = n_bins_per_block,
			breaks = x@breaks,
			centers = x@centers,
			positions = seq(x@bin_size, block_size, by = x@bin_size) - (block_size / 2)
		)
	}
) # compute_deviations
