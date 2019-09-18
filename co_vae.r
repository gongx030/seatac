
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
#gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 15
#gs <- 'Etv2_MEF_reprogramming'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 5
#gs <- 'D1_Dox_Etv2_all_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 5
#gs <- 'D1_Etv2_on_MEF_D1_D7Flk1pos_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 5
#gs <- 'MEF_NoDox+MEF_Dox_D7_Flk1pos+EB_Dox_D25_Flk1pos'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 10
#gs <- 'D1_Etv2_on_MEF_D7Flk1pos_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 5
gs <- 'D1_Etv2_on_MEF_D7Flk1pos_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 5

# --- prepare windows
source('analysis/seatac/helper.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train the model
latent_dim <- 2; epochs <- 20; type <- 'cvae'
#source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
#source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, min_reads_per_window = mr, type = type)
#source('analysis/seatac/helper.r'); save_model(model, model_dir)

# --- load the windows and the trained model
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], epochs)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)


# --- make the prediction
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)

R <- matrix(windows$num_reads, length(windows) / metadata(windows)$num_samples, metadata(windows)$num_samples)
i <- rowSums(R >= 10) == metadata(windows)$num_samples
i <- rep(i, metadata(windows)$num_samples)
gr <- windows[i]

library(irlba); u <- svd(mcols(gr)$latent)$u
k <- 2
set.seed(1); cls <- kmeans(u, k)$cluster
mcols(gr)$cluster <- cls

yy <- c(100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
for (h in 1:k){
	j <- mcols(gr)$cluster == h & gr$group == 1
	bg_counts <- colMeans(mcols(gr)$counts[gr$group == 1, ])
	xx <- as.matrix(colMeans(mcols(gr)$counts[j, ]) - bg_counts)
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('cluster=%d; n=%d', h, sum(j)),  breaks = c(-10, seq(quantile(xx, 0.01),quantile(xx, 0.99), length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)
	plot(colMeans(gr$coverage[j, ]))
	plot(colMeans(X[j, ]))
}

S <- matrix(gr$cluster, length(gr) / metadata(gr)$num_samples, metadata(gr)$num_samples)

# -----------------------------------------------------------------------------------
# [2019-08-27] Look at the V plot of Etv2 ChIP_seq peaks  at MEF and D7
# -----------------------------------------------------------------------------------
gs <- 'D1_Dox_Etv2_all_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; mr <- 15
source('analysis/seatac/windows.r'); peaks <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])


brg1_files <- c(
	'Alver' = '/panfs/roc/scratch/gongx030/datasets/dataset=Alver_version=20190407a/Brg1_summits.bed',
	'Chronis' = '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20190405a/Brg1_summits.bed',
	'NoDox' = '/panfs/roc/scratch/gongx030/datasets/dataset=Brg1_version=20190820a/MEF_NoDox_Brg1_summits.bed',
	'Dox_D1' = '/panfs/roc/scratch/gongx030/datasets/dataset=Brg1_version=20190820a/MEF_Dox_D1_Brg1_summits.bed'
)
devtools::load_all('packages/compbio'); brg1 <- lapply(brg1_files, function(brg1_file) macs2.read_summits(brg1_file))
mcols(windows)$brg1 <- do.call('cbind', lapply(brg1, function(x) resize(windows, fix = 'center', width = 200) %over% resize(x, fix = 'center', width = 200)))
colnames(mcols(windows)$brg1) <- names(brg1_files)

black <- !mcols(windows)$brg1[, 'Alver'] & !mcols(windows)$brg1[, 'Dox_D1']
red <- !mcols(windows)$brg1[, 'Alver'] & mcols(windows)$brg1[, 'Dox_D1']
green <- mcols(windows)$brg1[, 'Alver'] & !mcols(windows)$brg1[, 'Dox_D1']
blue <- mcols(windows)$brg1[, 'Alver'] & mcols(windows)$brg1[, 'Dox_D1']

yy <- c(100, 180, 247, 315, 473)

gr <- peaks

g <- 1; j <- mcols(gr)$group == g & mcols(gr)$num_reads >= 2 & red

bg_counts <- colMeans(mcols(gr)$counts[mcols(gr)$group == g, ])
xx <- as.matrix(colMeans(mcols(gr)$counts[j, ]) - bg_counts)
image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, breaks = c(-10, seq(quantile(xx, 0.01),quantile(xx, 0.99), length.out = 99), 10))
axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)


# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------
R <- matrix(windows$num_reads, length(windows) / metadata(windows)$num_samples, metadata(windows)$num_samples)
i <- rowSums(R >= 5) == metadata(windows)$num_samples
i <- rep(i, metadata(windows)$num_samples)
gr <- windows[i]

mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20170519a/MNase_treat_pileup.bw'
cvg <- rtracklayer::import(mnase_file, which = trim(reduce(gr)), as = 'RleList')
X <- as(as(cvg[gr], 'RleViews'), 'matrix')

library(irlba); u <- svd(mcols(gr)$latent)$u
k <- 2
set.seed(1); cls <- kmeans(u, k)$cluster
mcols(gr)$cluster <- cls

par(mfcol = c(3, k))
yy <- c(0, 100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
for (h in 1:k){
	j <- mcols(gr)$cluster == h & gr$group == 1
	bg_counts <- colMeans(mcols(gr)$counts[gr$group == 1, ])
	xx <- as.matrix(colMeans(mcols(gr)$counts[j, ]) - bg_counts)
	image(matrix(xx , metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('cluster=%d; n=%d', h, sum(j)),  breaks = c(-10, seq(quantile(xx, 0.01),quantile(xx, 0.99), length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2, lwd = 2)
	xx <- seq(0, 1, length.out = metadata(gr)$window_size)
	plot(xx, colMeans(gr$coverage[j, ]), xaxs="i", type = 'l', lwd = 4, xaxt = 'n', ylim = c(4, 5))
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'black', lty = 2, lwd = 2)
	plot(xx, colMeans(X[j, ]), xaxs="i", type = 'l', lwd = 4, xaxt = 'n', ylim = c(4, 5))
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'black', lty = 2, lwd = 2)
}

S <- matrix(gr$cluster, length(gr) / metadata(gr)$num_samples, metadata(gr)$num_samples)
gg <- sprintf('%s%s', S[, 1], S[, 2])

# --- Etv2 intensity

gr2 <- gr[gr$group == 1]
bw_file <-   '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_FE.bw'
bw_file <-   '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_NoDox_H3K27ac_FE.bw'
cvg <- rtracklayer::import(bw_file, which = trim(reduce(gr2)), as = 'RleList')
X <- as(as(cvg[gr2], 'RleViews'), 'matrix')

xx <- seq(0, 1, length.out = metadata(gr2)$window_size)
plot(xx, colMeans(X[gg == '12', ]), xaxs="i", type = 'l', lwd = 4, xaxt = 'n', ylim = c(0, 10), col = 'black')
lines(xx, colMeans(X[gg == '11', ]), lwd = 4, col = 'gray')
lines(xx, colMeans(X[gg == '21', ]), lwd = 4, col = 'red')
lines(xx, colMeans(X[gg == '22', ]), lwd = 4, col = 'gold')
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'black', lty = 2, lwd = 2)

bg_counts <- colMeans(mcols(gr2)$counts)
xx <- as.matrix(colMeans(mcols(gr2)$counts[gg == '22', ]) - bg_counts)
image(matrix(xx , metadata(gr2)$n_bins_per_window, metadata(gr2)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, breaks = c(-10, seq(quantile(xx, 0.01),quantile(xx, 0.99), length.out = 99), 10))





nuc2nfr <- gg == '12'
nfr2nuc <- gg == '21'
devtools::load_all('packages/compbio'); annot <- lapply(list(nuc2nfr, nfr2nuc), function(i) homer_annotatePeaks(gr[gr$group == 1][i], 'mm10'))







