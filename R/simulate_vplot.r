#' Simulate the V-plot
#' 
setMethod(
	'simulate_vplot',
	signature(
		x = 'missing'
	),
	function(x, ...){
		simulate_vplot_de_novo(...)
	}
)



#' Simulate the V-plot
#' 
setMethod(
	'simulate_vplot',
	signature(
		x = 'Vplots'
	),
	function(x, ...){

		param <- list(...)

		if (is.null(param$position)){

			z <- colSums(x$counts)
			z <- matrix(z, x@n_bins_per_window, x@n_intervals)
			browser()

		}
	}
)


simulate_vplot_de_novo <- function(
	n = 100, 
	mean_reads = 20,
	size_reads = 1.2,
	prob_nfr = 0.9,

	mean_distance_nfr = 0,
	sd_distance_nfr = 10,
	mean_distance_mono = 80,
	sd_distance_mono = 10,

	window_size = 640, 
	bin_size = 10, 
	fragment_size_range = c(0, 320), 
	fragment_size_interval = 10,
	mean_nfr_fragment_size = 100,
	mean_mono_fragment_size = (180 + 247) / 2,
	sd_mono_fragment_size = 30
){

	flog.info(sprintf('sampling %d counts from nbinom(mu=%.3f, size=.3f)', n, mean_reads, size_reads))
	n_reads <- rnbinom(n, mu = mean_reads, size = size_reads)
	ids <- rep(1:n, n_reads)
	m <- sum(n_reads)

	flog.info(sprintf('sampling state (nfr or mono_nucleosome) for each read, where prob of NFR is %.3f', m, prob_nfr))
	states <- sample(c('nfr', 'mono_nucleosome'), m, replace = TRUE, prob = c(prob_nfr, 1 - prob_nfr))  %>%
		factor()

	is_nfr <- states == 'nfr'
	is_mono <- states == 'mono_nucleosome'

	fs <- rep(0, m)
	flog.info(sprintf('sampling fragment size of %d NFR reads with exp(rate = %.3f)', sum(is_nfr), 1 / mean_nfr_fragment_size))
	fs[is_nfr] <- rexp(sum(is_nfr), rate = 1 / mean_nfr_fragment_size) %>%
		round()

	flog.info(sprintf('sampling fragment size of %d mono nucleosome reads with norm(mean = %.3f, sd = %.3f)', sum(is_mono), mean_mono_fragment_size, sd_mono_fragment_size))
	fs[is_mono] <- rnorm(sum(is_mono), mean = mean_mono_fragment_size, sd = sd_mono_fragment_size)

	breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
	centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
	fs <- as.numeric(cut(fs, breaks, include.lowest = TRUE)) # discretize the fragment size

	d <- rep(0, m)

	flog.info(sprintf('sampling distance-to-center of %d NFR reads with norm(mean=%.3f, sd=%.3f)', sum(is_nfr), mean_distance_nfr, sd_distance_nfr))
	d[is_nfr] <- rnorm(sum(is_nfr), mean = mean_distance_nfr, sd = sd_distance_nfr)

	flog.info(sprintf('sampling distance-to-center of %d mono nucleosome reads with norm(mean=%.3f, sd=%.3f)', sum(is_mono), mean_distance_mono , sd_distance_mono))
	d[is_mono] <- rnorm(sum(is_mono), mean = mean_distance_mono, sd = sd_distance_mono)

	d <- d * sample(c(1, -1), m, replace = TRUE)

	br <- seq(0, window_size, by = bin_size)
	d <- as.numeric(cut(d + window_size / 2, br, include.lowest = TRUE))


	n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval
	n_bins_per_window <- window_size / bin_size

	valid <- !is.na(fs) & !is.na(d)
	ids <- factor(ids[valid], 1:n)
	d <- factor(d[valid], 1:n_bins_per_window)
	fs <- factor(fs[valid], 1:n_intervals)

	counts <- table(ids, d, fs)
	counts <- as(counts, 'array')
	dim(counts) <- c(n, n_bins_per_window * n_intervals)
	counts <- as(counts, 'dgCMatrix')

	gr <- GRanges(seqnames = rep('chr1', n), range = IRanges(1, window_size))
	gr$counts <- counts
	new(
		'Vplots',
		gr,
		fragment_size_range  = as.integer(fragment_size_range),
		fragment_size_interval = as.integer(fragment_size_interval),
		bin_size = as.integer(bin_size),
		window_size = as.integer(window_size),
		n_intervals = as.integer(n_intervals),
		n_bins_per_window = as.integer(n_bins_per_window ),
		breaks = breaks,
		centers = centers,
		positions = seq(bin_size, window_size, by = bin_size) - (window_size / 2)
	)

}
