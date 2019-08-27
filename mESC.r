
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
gs <- 'Maza_mESC'; window_size <- 256; bin_size <- 2; fs <- c(50, 370, 5); genome <- 'mm10'

# --- prepare windows
source('analysis/seatac/helper.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train the model
latent_dim <- 2; epochs <- 20
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, min_reads_per_window = 10, min_reads_coverage = 10)
source('analysis/seatac/helper.r'); save_model(model, model_dir)

# --- load the windows and make the prediction
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size)
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, epochs = epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)

train <- mcols(windows)$num_reads >= 10 & mcols(windows)$mean_coverage >= 10 
gr <- windows[train]
library(irlba); u <- svd(mcols(gr)$latent, nu = 2, nv = 1)$u
set.seed(1); cls <- kmeans(u, 4)$cluster
mcols(gr)$cluster <- cls

par(mfcol = c(3, 4))
k <- max(cls)
yy <- c(50, 200, 400, 600, 670)
lapply(1:k, function(h){
	i <- mcols(gr)$cluster == h
	image(matrix(as.matrix(colMeans(mcols(gr)$counts[i, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('C%d(%d)', h, sum(i)))
	axis(2, at = (yy - 50) / (670 - 50), label = yy)
	abline(v = c(0.5 - 180/metadata(gr)$window_size , 0.5, 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	axis(1, at = c(0, 0.5, 1))
	plot(colMeans(mcols(gr)$coverage[i,]), type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'coverage')
	x <- 1:metadata(gr)$window_size
	y <- colMeans(mcols(gr)$nucleosome_score[i, ])
	plot(x, y, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'NCP score')
	lo <- loess(y~x)
	lines(predict(lo), col='red', lwd=2)
	abline(v = metadata(gr)$window_size / 2)
})





