
library(tensorflow)
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)

# -----------------------------------------------------------------------------------
# [2019-08-27] mESC data
# ------------------------------------------------------------------------------------
gs <- 'Maza_mESC'; window_size <- 640; bin_size <- 10; fs <- c(50, 690, 10); genome <- 'mm10'
gs <- 'Maza_mESC'; window_size <- 480; bin_size <- 5; fs <- c(50, 690, 10); genome <- 'mm10'
gs <- 'Maza_mESC'; window_size <- 480; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
gs <- 'Maza_mESC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'

# --- prepare windows
source('analysis/seatac/helper.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train the model
latent_dim <- 2; epochs <- 20
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, min_reads_per_window = 30)
source('analysis/seatac/helper.r'); save_model(model, model_dir)

# --- load the windows and make the prediction
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size)
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, epochs = epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)

train <- mcols(windows)$num_reads >= 30 & mcols(windows)$mean_coverage >= 0 
gr <- windows[train]
library(irlba); u <- svd(mcols(gr)$latent, nu = 3, nv = 1)$u
set.seed(1); cls <- kmeans(u, 4)$cluster
#set.seed(1); cls <- kmeans(mcols(gr)$latent, 5)$cluster
mcols(gr)$cluster <- cls

par(mfcol = c(4, 5))
k <- max(cls)
yy <- c(100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
lapply(1:k, function(h){
	i <- mcols(gr)$cluster == h
#	xx <- as.matrix(colMeans(mcols(gr)$counts[i, ])) - colMeans(mcols(gr)$counts)
	xx <- as.matrix(colMeans(mcols(gr)$counts[i, ]))
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('C%d(%d)', h, sum(i)))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - 180/metadata(gr)$window_size , 0.5, 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	axis(1, at = c(0, 0.5, 1))
	x <- mcols(gr)$counts[i, ] %>% as.matrix()
	dim(x) <- c(nrow(x), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	nfr <- rowMeans(colSums(x[, , which(breaks <= 100)], dim = 1))
	mono_nucleosome <- rowMeans(colSums(x[, , which(breaks >= 180 & breaks <= 247)], dim = 1))
#	plot(mono_nucleosome / nfr, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'Mono nucleosome')
	plot(mono_nucleosome , type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'Mono nucleosome')
	plot(colMeans(mcols(gr)$coverage[i,]), type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'coverage')
	x <- 1:metadata(gr)$window_size
	y <- colMeans(mcols(gr)$nucleosome_score[i, ])
	plot(x, y, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'NCP score')
	lo <- loess(y~x, span = 0.1)
	lines(predict(lo), col='red', lwd=2)
	abline(v = metadata(gr)$window_size / 2)
})



x <- mcols(gr)$counts %>% as.matrix()
dim(x) <- c(nrow(x), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
x_mean <- colSums(x, dim = 1)


