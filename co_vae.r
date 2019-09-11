
library(tensorflow)
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)

# -----------------------------------------------------------------------------------
# [2019-08-27] combined mESC and MEF data
# -----------------------------------------------------------------------------------
#gs <- 'MEF_ESC'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 30
#gs <- 'MEF_ESC_with_MEF_Peaks'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 2
#gs <- 'MEF_ESC_with_MEF_Peaks'; window_size <- 480; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 2
#gs <- 'MEF_ESC_with_MEF_Peaks'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 10); genome <- 'mm10'; mr <- 2
#gs <- 'MEF_ESC_with_MEF_Peaks'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 2
#gs <- 'D1_Dox_Etv2_on_D0+D1_MEF'; window_size <- 480; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 2
#gs <- 'D2_Dox_Etv2_on_D0+D2_MEF'; window_size <- 480; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 2
gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 5

# --- prepare windows
source('analysis/seatac/helper.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train the model
#latent_dim <- 2; epochs <- 20
latent_dim <- 10; epochs <- 50
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
devtools::load_all('analysis/seatac/packages/seatac'); windows <- filter_windows(windows, min_reads_per_window = mr)
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs)
source('analysis/seatac/helper.r'); save_model(model, model_dir)

# --- load the windows and the trained model
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)


# --- make the prediction
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)
gr <- windows

s <- mcols(gr)$group == 1
library(irlba); u <- svd(mcols(gr)$latent[s, ])$u
k <- 5
set.seed(1); cls <- kmeans(u, k)$cluster
mcols(gr)$cluster <- cls

par(mfrow = c(metadata(windows)$num_samples, k))
yy <- c(100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
for (i in 1:metadata(windows)$num_samples){
	for (h in 1:k){
		j <- mcols(gr)$cluster == h & mcols(gr)$group == i
		bg_counts <- colMeans(mcols(gr)$counts[mcols(gr)$group == i, ])
#		xx <- as.matrix(colMeans(mcols(gr)$counts[j, ]) - bg_counts)
		xx <- as.matrix(colMeans(mcols(gr)$counts[j, ]))
		image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('group=%d; cluster=%d; n=%d', i, h, sum(j)),  breaks = c(-10, seq(quantile(xx, 0.01),quantile(xx, 0.99), length.out = 99), 10))
		axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
		abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	}
}

