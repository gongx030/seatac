
library(tensorflow)
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)

# -----------------------------------------------------------------------------------
# [2019-08-27] mESC data
# ------------------------------------------------------------------------------------
#gs <- 'Maza_mESC'; window_size <- 640; bin_size <- 10; fs <- c(50, 690, 10); genome <- 'mm10'; mr <- 30
#gs <- 'Maza_mESC'; window_size <- 480; bin_size <- 5; fs <- c(50, 690, 10); genome <- 'mm10'; mr <- 30
#gs <- 'Maza_mESC'; window_size <- 480; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 30
#gs <- 'Maza_mESC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 30
gs <- 'Maza_mESC'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 30
#gs <- 'Maza_mESC'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 30
#gs <- 'Maza_mESC'; window_size <- 200; bin_size <- 1; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 10

# --- prepare windows
source('analysis/seatac/helper.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train the model
latent_dim <- 2; epochs <- 20
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, min_reads_per_window = mr)
source('analysis/seatac/helper.r'); save_model(model, model_dir)

# --- load the windows and the trained model
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)

# --- make the prediction
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)
train <- mcols(windows)$num_reads >= mr & mcols(windows)$mean_coverage >= 5 
gr <- windows[train]

library(irlba); u <- svd(mcols(gr)$latent)$u
set.seed(1); cls <- kmeans(u, 4)$cluster
mcols(gr)$cluster <- cls

par(mfcol = c(6, 7))
k <- max(cls)
yy <- c(100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
lapply(1:k, function(h){
	i <- mcols(gr)$cluster == h
	xx <- as.matrix(colMeans(mcols(gr)$counts[i, ])) - colMeans(mcols(gr)$counts)
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('C%d(%d)', h, sum(i)), breaks = c(-10, seq(-0.01, 0.01, length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - 180/metadata(gr)$window_size , 0.5, 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	axis(1, at = c(0, 0.5, 1))

	x <- as.matrix(t(t(mcols(gr)$counts[i, ]) - colMeans(mcols(gr)$counts)))
	dim(x) <- c(nrow(x), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	nfr <- rowMeans(colSums(x[, , which(breaks <= 100)], dim = 1))
	mono_nucleosome <- rowMeans(colSums(x[, , which(breaks >= 180 & breaks <= 247)], dim = 1))
	plot(mono_nucleosome - nfr, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'Mono nucleosome')

	xx <- as.matrix(colMeans(mcols(gr)$counts[i, ]))
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('C%d(%d)', h, sum(i)), breaks = c(-1, seq(0, 0.05, length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - 180/metadata(gr)$window_size , 0.5, 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	axis(1, at = c(0, 0.5, 1))

	plot(colMeans(mcols(gr)$coverage[i,]), type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'coverage')
	plot(colMeans(mcols(gr)$nucleoatac_signal[i,]), type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'Nucleoatac')

	x <- 1:metadata(gr)$window_size

	y <- colMeans(mcols(gr)$nucleosome_score[i, ])
	plot(x, y, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'NCP score')
	lo <- loess(y~x, span = 0.1)
	lines(predict(lo), col='red', lwd=2)
	abline(v = metadata(gr)$window_size / 2)

})


# -----------------------------------------------------------------------------------
# [2019-09-04] Testing the trained model on MEF
# ------------------------------------------------------------------------------------
gs <- 'Maza_mESC'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 30
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)

gs <- 'MEF_NoDox'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 30
source('analysis/seatac/helper.r'); win_test <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
devtools::load_all('analysis/seatac/packages/seatac'); win_test <- model %>% predict(win_test)
train <- mcols(win_test)$num_reads >= mr & mcols(win_test)$mean_coverage >= 5 
gr <- win_test[train]

library(irlba); u <- svd(mcols(gr)$latent)$u
set.seed(1); cls <- kmeans(u, 4)$cluster
mcols(gr)$cluster <- cls

par(mfcol = c(6, 7))
k <- max(cls)
yy <- c(100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
lapply(1:k, function(h){

	i <- mcols(gr)$cluster == h
	bg_counts <- colMeans(mcols(gr)$counts)
	xx <- as.matrix(colMeans(mcols(gr)$counts[i, ]) - bg_counts)
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('C%d(%d)', h, sum(i)), breaks = c(-10, seq(-0.05, 0.05, length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - 180/metadata(gr)$window_size , 0.5, 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	axis(1, at = c(0, 0.5, 1))

	x <- as.matrix(t(t(mcols(gr)$counts[i, ]) - bg_counts))
	dim(x) <- c(nrow(x), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	nfr <- rowMeans(colSums(x[, , which(breaks <= 100)], dim = 1))
	mono_nucleosome <- rowMeans(colSums(x[, , which(breaks >= 180 & breaks <= 247)], dim = 1))
	plot(mono_nucleosome - nfr, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'Mono nucleosome')

	xx <- as.matrix(colMeans(mcols(gr)$counts[i, ]))
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('C%d(%d)', h, sum(i)), breaks = c(-1, seq(0, 0.05, length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - 180/metadata(gr)$window_size , 0.5, 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	axis(1, at = c(0, 0.5, 1))

	plot(colMeans(mcols(gr)$coverage[i,]), type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'coverage')

	x <- 1:metadata(gr)$window_size

	y <- colMeans(mcols(gr)$nucleosome_score[i, ])
	plot(x, y, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'Nucleosome score')
	lo <- loess(y~x, span = 0.1)
	lines(predict(lo), col='red', lwd=2)
	abline(v = metadata(gr)$window_size / 2)

	y3 <- colMeans(mcols(gr)$fitted_weighted_nucleosome_score[i, ])
	plot(x, y3, type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'fitted weighted NCP score')
	lo <- loess(y3~x, span = 0.1)
	lines(predict(lo), col='red', lwd=2)
})


