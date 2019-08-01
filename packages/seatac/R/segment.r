segment <- function(x, k = 2, p0 = 0.9){

	num_steps <- metadata(x)$n_bins_per_window
  feature_dim <- ncol(mcols(x)$fitted_counts) / num_steps

	breaks <- seq(metadata(x)$fragment_size_range[1], metadata(x)$fragment_size_range[2], by = metadata(x)$fragment_size_interval)

	is_nfr <- breaks > 0 & breaks <= 100
	is_mono_nucleosome <- breaks > 180 & breaks <= 247	
	is_di_nucleosome <- breaks > 315 & breaks <= 473
	is_tri_nucleosome <- breaks > 558 & breaks <= 615

	X <- mcols(x)$fitted_counts %>%
		tf$reshape(shape(length(x), feature_dim, num_steps))
		
	X_sum <- tf$reduce_sum(X, axis = 1L, keepdims = TRUE)
	X <- tf$div(X, X_sum) 
	X <- X %>% 
		tf$transpose(perm = c(0L, 2L, 1L)) %>% 
		tf$reshape(shape(length(x) * num_steps, feature_dim)) %>%
		as.matrix()

	ntimes <- rep(num_steps, length(x))

	y <- cbind(
		nfr = rowSums(X[, is_nfr]),
		mono_nucleosome = rowSums(X[, is_mono_nucleosome]),
		di_nucleosome = rowSums(X[, is_di_nucleosome]),
		tri_nucleosome = rowSums(X[, is_tri_nucleosome])
	)

	trstart <- matrix((1 - p0)/ (k - 1), k, k)
	diag(trstart) <- p0
	transition <- list()

	rModels <- list()
	for (i in 1:k){
		rModels[[i]] <- list(MVNresponse(y ~ 1))
		transition[[i]] <- transInit(~1, nstates = k, data = data.frame(1), pstart = trstart[i, ])
	}

	instart <- rep(1 / k, k)

	inMod <- transInit(~1, ns = k, ps = instart, data = data.frame(rep(1, length(x))))

	mod <- makeDepmix(response = rModels, transition = transition, prior = inMod, ntimes = ntimes)
	fm <- depmixS4::fit(mod)

	state <- posterior(fm)[, 'state']

	h <- order(sapply(1:k, function(i) mean(y[state == i, 'nfr'])))
	state <- h[state]
	
	mcols(x)$state <- matrix(state, length(x), num_steps, byrow = TRUE)
	x

} # segment


call_nucleosome <- function(x){
	length(x)
	browser()
	lapply(1:length(x), function(i){
		rle(mcols(x)$state[i, ])
	})
} # call_nucleosome
