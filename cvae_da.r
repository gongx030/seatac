
library(tensorflow)
tfe_enable_eager_execution(device_policy = 'silent')
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)


#gs <- 'MEF_NoDox'; flanking <- 1000; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; genome2 <- BSgenome.Mmusculus.UCSC.mm10
gs_train <- 'MEF_NoDox_Flk1pos_D7'; flanking <- 1000; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; genome2 <- BSgenome.Mmusculus.UCSC.mm10

# --- prepare windows
source('analysis/seatac/prepare_windows.r'); windows <- prepare_windows(gs_train, flanking, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs_train, flanking, bin_size, fs[1:2], fs[3])

# --- load training windows
source('analysis/seatac/windows.r'); windows <- load_windows(gs_train, flanking, bin_size, fs[1:2], fs[3])

# --- train model
#latent_dim <- 5; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 50; bs <- 64; spe <- 20
#latent_dim <- 10; sequence_dim <- 32; n_components <- 15; mr <- 50; window_size <- 640; epochs <- 50; bs <- 64; spe <- 20
latent_dim <- 10; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 5; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 10; bs <- 64; spe <- 20
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs_train, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, n_components, mr, flanking, epochs, bs, spe)

devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, window_size = window_size, min_reads_per_window = mr, epochs = epochs, batch_size = bs, steps_per_epoch = spe, sequence_dim = sequence_dim, n_components = n_components, genome = genome2)
devtools::load_all('analysis/seatac/packages/seatac'); save_model(model, model_dir)

# --- load models
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)

# --- prepare testing windows
gs_test <- 'motif_on_MEF_NoDox_Flk1pos_D7'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; genome2 <- BSgenome.Mmusculus.UCSC.mm10
#gs_test <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'; genome2 <- BSgenome.Mmusculus.UCSC.mm10

source('analysis/seatac/prepare_windows.r'); windows <- prepare_windows(gs_test, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs_test, window_size, bin_size, fs[1:2], fs[3])

# --- load testing windows
source('analysis/seatac/windows.r'); windows <- load_windows(gs_test, window_size, bin_size, fs[1:2], fs[3])

# --- make prediction
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows[windows$group == 1], batch_size = 128)

core_width <- 180
bin_size <- metadata(model$data)$bin_size
core_bin <- round(metadata(model$data)$n_bins_per_window / 2 - (core_width / 2) / bin_size):round(metadata(model$data)$n_bins_per_window / 2 + (core_width / 2) / bin_size)
nfr <- rowMeans(model$data$nfr[, core_bin])
mono <- rowMeans(model$data$mono_nucleosome[, core_bin])
r <- log(mono + 1) - log(nfr + 1)

r2 <- log(rowMeans(gr$mono_nucleosome[, core_bin]) + 1) - log(rowMeans(gr$nfr[, core_bin]) + 1)
z_score <- (r2 - mean(r)) / sd(r)

k <- 3
cls <- rep(3, length(gr))
cls[z_score > qnorm(0.6)] <- 1
cls[z_score < qnorm(0.4)] <- 2
mcols(gr)$cluster <- cls

#Z <- gr$latent
Z <- gr$latent +  gr$h
#u <- prcomp(Z)$x
#set.seed(1); cls <- kmeans(Z, k)$cluster; mcols(gr)$cluster <- cls

par(mfcol = c(3, k + 1), mar = c(2, 5, 2, 1))
yy <- c(50, 100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
for (h in 1:k){

	j <- mcols(gr)$cluster == h 
	xx <- as.matrix(colMeans(mcols(gr)$counts[j, ]))
	bg_counts <- matrix(colMeans(mcols(gr)$counts), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	X <- matrix(xx, metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	X <- X - bg_counts
	image(X, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = sprintf('cluster=%d; n=%d', h, sum(j)),  breaks = c(-10, seq(quantile(X, 0.01),quantile(X, 0.99), length.out = 99), 10))
	axis(2, at = (yy - fs[1]) / (fs[2] - fs[1]), label = yy)
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'yellow', lty = 2)

	xx <- seq(0, 1, length.out = metadata(gr)$n_bins_per_window)
	d <- data.frame(x = xx, y = log(rowMeans(X[, metadata(gr)$mono_nucleosome]) + 0.1) - log(rowMeans(X[, metadata(gr)$nfr]) + 0.1))
	m <- loess(y ~ x, data = d, span = 0.4)
	plot(xx, predict(m), lwd = 3, col = 'blue', type = 'l', xaxt = 'n', xaxs = 'i', ylim = c(-0.15, 0.15), xlab = '', ylab = 'Likelihood ratio')
	abline(h = 0, lwd = 1, lty = 2)
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)

	xx <- seq(0, 1, length.out = metadata(gr)$window_size)
#	plot(xx, colMeans(gr$nucleosome_score[j, ]), xaxt = 'n', xaxs = 'i', ylim = range(colMeans(gr$nucleosome_score)) * c(0.9, 1.1), ylab = 'Read density', xlab = '')
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)

#	plot(xx, colMeans(gr$nucleoatac_signal[j, ]), xaxt = 'n', xaxs = 'i', ylab = 'Read density', xlab = '', ylim = c(0, 0.3))
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
}

plot(Z[, 1:2], pch = 21, bg = cls)

pairs(Z, bg = 21, col = cls + 1)


library(Rtsne); set.seed(1); y_tsne <- Rtsne(u, check_duplicates = FALSE)$Y
plot(y_tsne, col = cls, bg = 21)


# ----------------------------------------------------------------------------------------
# [2019-11-01] Predicting the nucleosome/nfr conversion between MEF and D7 Flk1+ cells
# for every motif
# ----------------------------------------------------------------------------------------
gs_train <- 'MEF_NoDox_Flk1pos_D7'; genome <- 'mm10'; latent_dim <- 10; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); flanking <- 1000; epochs <- 100; bs <- 64; spe <- 20
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs_train, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, n_components, mr, flanking, epochs, bs, spe)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)

# --- prepare all motifs and their positions
library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(motifmatchr)
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')
se <- SummarizedExperiment(rowRanges = model$input_data)
peaks <- matchMotifs(get(motif.set), se, genome = 'mm10', out = 'positions')
motifs <- rep(names(peaks), sapply(peaks, length))
peaks <- Reduce('c', peaks)
peaks$motifs <- motifs
peaks <- resize(peaks, fix = 'center', width = window_size)
peaks <- peaks[order(peaks)]
peaks <- unique(peaks)
devtools::load_all('packages/compbio'); peaks <- add.seqinfo(peaks, 'mm10')
se_file <- sprintf('%s/MEF_NoDox_Flk1pos_D7_ALL_motifs.rds', project_dir('seatac'))
saveRDS(peaks, se_file)



  ga <- read_bam(bam_files, peaks, genome = genome2, expand = window_size * 2)
  peaks <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size, fragment_size_range = fragment_size_range, fragment_size_interval = fragment_size_interval)


for (i in 2:length(get(motif.set))){
	pwm <- get(motif.set)[[i]]	# a PWMatrix object
	window_file <- sprintf('%s/motif_windows/gs_train=%s_window_size=%d_motif=%s.rds', project_dir('seatac'), gs_train, window_size, pwm@name)
	flog.info(sprintf('%d/%d writing %s', i, length(get(motif.set)), window_file))
	saveRDS(windows, window_file)
}


se <- SummarizedExperiment(rowRanges = peaks)
motif_ix <- matchMotifs(get(motif.set)[hs], se, genome = 'mm10', out = 'positions')


