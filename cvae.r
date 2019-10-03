

library(tensorflow)
tfe_enable_eager_execution(device_policy = 'silent')
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)

# -----------------------------------------------------------------------------------
# [2019-08-27] combined mESC and MEF data
# -----------------------------------------------------------------------------------
#gs <- 'MEF_NoDox'; window_size <- 320; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'MEF_NoDox'; window_size <- 64; bin_size <- 2; fs <- c(50, 370, 5); genome <- 'mm10'
gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 320; bin_size <- 10; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'Maza_mESC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'Maza_mESC'; window_size <- 320; bin_size <- 2; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'GM12878'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'hg19'
#gs <- 'D1_Dox_Etv2'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'

# --- prepare windows
source('analysis/seatac/prepare_windows.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train model
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
latent_dim <- 2; epochs <- 25; sequence_dim <- 32; type <- 'gmm_cvae'; n_components <- 15
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, sequence_dim = sequence_dim, type = type, n_components = n_components)
#source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim)
#devtools::load_all('analysis/seatac/packages/seatac'); save_model(model, model_dir)

# --- load model and make prediction
#source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim)
#devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows)

Z <- gr$latent
k <- 5
u <- svd(Z)$u
set.seed(1); cls <- kmeans(Z, k)$cluster
mcols(gr)$cluster <- cls


par(mfcol = c(4, k))
yy <- c(50, 100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
for (h in 1:k){
	j <- mcols(gr)$cluster == h 
	bg_counts <- colMeans(mcols(gr)$counts); x <- as.matrix(colMeans(mcols(gr)$counts[j, ]) - bg_counts)
	X <- matrix(x, metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	image(X, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('cluster=%d; n=%d', h, sum(j)),  breaks = c(-10, seq(quantile(x, 0.01),quantile(x, 0.99), length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)

#	xx <- seq(0, 1, length.out = metadata(gr)$n_bins_per_window)

#	Y <- t(t(gr$mono_nucleosome[j, ]) - colMeans(gr$mono_nucleosome))
#	Y <- (Y - rowMins(Y)) / (rowMaxs(Y) - rowMins(Y))
#	Y2 <- t(t(gr$nfr[j, ]) - colMeans(gr$nfr))
#	Y2 <- (Y2 - rowMins(Y2)) / (rowMaxs(Y2) - rowMins(Y2))

#	xx <- seq(0, 1, length.out = metadata(gr)$n_bins_per_window)
#	d <- data.frame(x = xx, y = colMeans(log(Y + 0.1) - log(Y2 + 0.1)))
#	d <- data.frame(x = xx, y = log(rowMeans(X[, metadata(gr)$mono_nucleosome]) + 0.1) - log(rowMeans(X[, metadata(gr)$nfr]) + 0.1))
#	m <- loess(y ~ x, data = d, span = 0.2)
#	plot(xx, predict(m), lwd = 3, col = 'blue', type = 'l')
#	abline(h = 0, lwd = 1, lty = 2)
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)

	xx <- seq(0, 1, length.out = metadata(gr)$window_size)
	plot(xx, colMeans(gr$coverage[j, ]))
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
	plot(xx, colMeans(gr$nucleosome_score[j, ]))
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
	plot(xx, colMeans(gr$nucleoatac_signal[j, ]))
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
}

par(mfrow = c(4, 6)); xx <- lapply(1:24, function(i) {plot(gr$chipseq[, 1][i, ], ylim = c(0, 1)); lines(gr$predicted_chipseq[, 1][i, ], col = 'red')})


Y <- t(t(gr$mono_nucleosome) - colMeans(gr$mono_nucleosome))
Y <- (Y - rowMins(Y)) / (rowMaxs(Y) - rowMins(Y))
Y2 <- t(t(gr$nfr) - colMeans(gr$nfr))
Y2 <- (Y2 - rowMins(Y2)) / (rowMaxs(Y2) - rowMins(Y2))
R <- log(Y + 0.1) - log(Y2 + 0.1)

j <- order(rowSums(R[, c(64 - 10):(64 + 10)]), decreasing = FALSE)[1:5000]



devtools::load_all('analysis/seatac/packages/seatac'); bins <- model %>% scan_motifs(gr, genome = BSgenome.Mmusculus.UCSC.mm10)

k <- 7
b <- bins[bins$window_id %in% gr$window_id[gr$cluster == k]]

A <- sparseMatrix(i = 1:length(b),j = b$bin_id, dims = c(length(b), metadata(gr)$window_size - model$kernel_size + 1))
A <- t(b$motifs) %*% A

i <- b$bin_id == 150

ss <- AlignSeqs(b[b$motifs[, 42] > 0]$sequence, verbose = FALSE, processors = 8)
pfm <- consensusMatrix(ss)

pfm <- pfm[c('A', 'C', 'G', 'T'), ]
pfm <- pfm + 5
icm <- pfm %*% diag(1 / colSums(pfm))
seqLogo::seqLogo(icm)

ss <- AlignSeqs(subseq(gr[gr$cluster == 7]$sequence, 200 - 15, 200 + 15), verbose = FALSE, processors = 8)
pfm <- consensusMatrix(ss)

pfm <- pfm[c('A', 'C', 'G', 'T'), ]
pfm <- pfm + 1
icm <- pfm %*% diag(1 / colSums(pfm))
seqLogo::seqLogo(icm)

