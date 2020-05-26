#' Build a logistic regression model for fragment sizes
#' 
#' Build a logistic regression model for classiying fragment sizes into NFR (nucleosome free region, below 100 bp) and
#' mono nucleosome (between 180 and 247 bp)
#'
#' @param x a GRange object
#' @param n number of sampling
#' 
#' @return a numeric vector
#'
build_fragment_size_model <- function(x, n = 1e5, bootstrap = 10){

	y <- colSums(x$smoothed_counts > 0) %>%
	 	matrix(
			metadata(x)$n_bins_per_window,
			metadata(x)$n_intervals,
		) %>%
	  colSums()

	d <- bplapply(1:bootstrap, function(i){

		di <- data.frame(
			fragment_size = sample(metadata(x)$centers, n, replace = TRUE, prob = y)
		)

		m <- flexmix(
			fragment_size ~ 1, 
			data = di, 
			k = 2,
			model = list(FLXMRglm(family = 'Gamma'), FLXMRglm(family = 'gaussian'))
		)

		cbind(di, cluster = clusters(m) - 1)

	})

	P <- do.call('cbind', bplapply(seq_len(bootstrap), function(i){
		model <- glm(cluster ~ fragment_size, d[[i]], family = binomial(link = 'logit'))
		p <- stats::predict(model, data.frame(fragment_size = metadata(x)$centers), type = 'response')
		if (p < 0.5)
			p <- 1 - p
		p
	}))

	rowSums(P) / bootstrap

} #
