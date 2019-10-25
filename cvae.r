

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
#gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 320; bin_size <- 10; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 10; fs <- c(50, 370, 10); genome <- 'mm10'
#gs <- 'Maza_mESC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'Maza_mESC'; window_size <- 320; bin_size <- 2; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'GM12878'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'hg19'
#gs <- 'D1_Dox_Etv2'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D1_Dox_Etv2_all_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'Etv2_MEF_reprogramming'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D1_Dox_Etv2_D1_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D2_Dox_Etv2_D2_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D7_Dox_Etv2_D7_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'D7_Dox_Etv2_D7Flk1pos_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'EB25_Etv2_Flk1pos_ATAC'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
#gs <- 'GM12878_STAT3'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'hg19'
#gs <- 'GM12878_ETV6'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'hg19'
#gs <- 'GM12878_STAT3'; window_size <- 640; bin_size <- 10; fs <- c(50, 370, 10); genome <- 'hg19'
gs <- 'PHF7_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'

# --- prepare windows
source('analysis/seatac/prepare_windows.r'); windows <- prepare_windows(gs, window_size, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs, window_size, bin_size, fs[1:2], fs[3])

# --- train model
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
latent_dim <- 5; epochs <- 50; sequence_dim <- 32; type <- 'gmm_cvae'; n_components <- 15; mr <- 20; bs <- 64
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, sequence_dim = sequence_dim, type = type, n_components = n_components, min_reads_per_window = mr, batch_size = bs)
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, type, n_components, mr)
devtools::load_all('analysis/seatac/packages/seatac'); save_model(model, model_dir)

# --- load model and make prediction
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, type, n_components, mr)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows, batch_size = 512)

#source('analysis/seatac/windows.r'); save_windows(gr, gs, window_size, bin_size, fs[1:2], fs[3])
#source('analysis/seatac/windows.r'); gr <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])


Z <- gr$latent
k <- 4
u <- prcomp(Z)$x
set.seed(1); cls <- kmeans(Z, k)$cluster; mcols(gr)$cluster <- cls


par(mfcol = c(3, k + 1))
yy <- c(50, 100, 180, 247, 315, 473)
breaks <- seq(fs[1], fs[2], by = fs[3])
for (h in 1:k){

	j <- mcols(gr)$cluster == h 
	x <- as.matrix(colMeans(mcols(gr)$counts[j, ]))
	bg_counts <- matrix(colMeans(mcols(gr)$counts), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
	X <- matrix(x, metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
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

#	plot(xx, colMeans(gr$seatac_nucleosome_ratio[j, ]), type = 'l', xaxt = 'n', xaxs = 'i', lwd = 3, col = 'darkblue', ylab = '')
#	abline(h = 0, lwd = 1, lty = 2)
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)

#	plot(xx, colMeans(gr$mono_nucleosome[j, ]) - colMeans(gr$mono_nucleosome), type = 'l', xaxt = 'n', xaxs = 'i', lwd = 3, col = 'darkblue', ylab = '')
#	abline(h = 0, lwd = 1, lty = 2)
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)

#	plot(xx, colMeans(gr$seatac_nucleosome_log10pvalue[j, ]), type = 'l', xaxt = 'n', xaxs = 'i', lwd = 3, col = 'darkblue', ylab = '')
#	abline(h = 0, lwd = 1, lty = 2)
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)

	xx <- seq(0, 1, length.out = metadata(gr)$window_size)
#	plot(xx, colMeans(gr$coverage[j, ]), col = h, xaxt = 'n', xaxs = 'i')
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
	plot(xx, colMeans(gr$nucleosome_score[j, ]), xaxt = 'n', xaxs = 'i', ylim = range(colMeans(gr$nucleosome_score)) * c(0.9, 1.1), ylab = 'Read density', xlab = '')
	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
#	plot(xx, colMeans(gr$nucleoatac_signal[j, ]), xaxt = 'n', xaxs = 'i', ylim = c(0.0, 0.5))
#	plot(xx, colMeans(gr$nucleoatac_signal[j, ]), xaxt = 'n', xaxs = 'i')
#	abline(v = c(0.5 - (1:4) * 180 /metadata(gr)$window_size , 0.5, (1:4) * 180/ metadata(gr)$window_size + 0.5), col = 'red', lty = 2)
}
plot(u[, 1:2], bg = 21, col = cls)


library(Rtsne); set.seed(1); y_tsne <- Rtsne(u, check_duplicates = FALSE)$Y
plot(y_tsne, col = cls, bg = 21)



# --- Look at the correlation between ground truth and inferred
bg_counts <- colMeans(mcols(gr)$counts); x <- as.matrix(colMeans(mcols(gr)$counts - bg_counts)
X <- matrix(x, metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)

par(mfcol = c(4, 6))
i <- 1
plot(gr$seatac_nucleosome_ratio[i, ])
plot(gr$seatac_nucleosome_log10pvalue[i, ])
plot(gr$nucleosome_score[i, ])
plot(gr$nucleoatac_signal[i, ])

A <- sparseMatrix(i = 1:metadata(gr)$window_size, j = rep(1:metadata(gr)$n_bins_per_window, each = metadata(gr)$bin_size), dims = c(metadata(gr)$window_size, metadata(gr)$n_bins_per_window))
A <- A %*% Diagonal(x = 1 / colSums(A))


xx <-  gr$seatac_nucleosome_log10pvalue
yy <- as.matrix(gr$nucleosome_score %*% A)
smoothScatter(xx[, 64], yy[, 64], xlim = c(0, 2), ylim = c(0, 5))
cor(xx[, 64], yy[, 64])

xx2 <- gr$nucleosome_score
yy2 <- gr$nucleoatac_signal
plot(xx2[, 320], yy2[, 320], xlim = c(0, 5))

smoothScatter(c(gr$seatac_nucleosome_log10pvalue), c(as.matrix(gr$nucleosome_score %*% A)), ylim = c(0, 5))
smoothScatter(c(gr$nucleosome_score), c(gr$nucleoatac_signal), xlim = c(0, 5), ylim = c(0, 5))


# -----------------------------------------------------------------------------------
# [2019-08-27] combined mESC and MEF data
# -----------------------------------------------------------------------------------
gs <- 'D1_Dox_Etv2_on_MEF'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
latent_dim <- 5; epochs <- 50; sequence_dim <- 32; type <- 'gmm_cvae'; n_components <- 15; mr <- 20
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, type, n_components, mr)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)


gs <- 'D1_Dox_Etv2_all_ATAC'
source('analysis/seatac/windows.r'); windows <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows)

gs <- 'D1_Dox_Etv2_on_MEF'
source('analysis/seatac/windows.r'); windows2 <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])
devtools::load_all('analysis/seatac/packages/seatac'); gr2 <- model %>% predict(windows2)

k <- 4
set.seed(1); cls2 <- kmeans(gr2$latent, k)$cluster

library(FNN)
cls <- knn(gr2$latent, gr$latent, factor(cls), k = 1)
A <- matrix(cls, nrow = length(cls) / 5, 5)
class(A) <- 'numeric'
mcols(gr)$cluster <- cls


# -----------------------------------------------------------------------------------
# [2019-10-07] Generate a universal GRange objet for all ENCODE ChIP-seq data for GM12878
# and K562, and use NucleoATAC to call NFR and nucleosome
# -----------------------------------------------------------------------------------
devtools::load_all('packages/compbio');
#dataset <- 'dataset=GM12878_version=20191007a'; genome <- 'hg19'
dataset <- 'dataset=K562_version=20191007a'; genome <- 'hg19'
d <- read_dataset(dataset)
grs <- lapply(1:nrow(d), function(i){
	x <- read.table(gzfile(d[i, 'bed_file']), header = FALSE, sep = '\t')
	GRanges(seqnames = x[, 1], range = IRanges(x[, 2], x[, 3]))
})
gr <- Reduce('c', (lapply(grs, function(gr) resize(gr, fix = 'center', width = 2000))))
gr <- reduce(gr)
gr <- add.seqinfo(gr, genome)
blacklist_file <- '/panfs/roc/scratch/gongx030/datasets/datasets=blacklists_version=20190827a/hg38.blacklist.bed.gz'
blacklist <- read.table(gzfile(blacklist_file), sep = '\t')
blacklist <- GRanges(seqnames = blacklist[, 1], range = IRanges(blacklist[, 2], blacklist[, 3]))
i <- gr %over% blacklist
gr <- gr[!i]
peak_file <- sprintf('%s/all_ChIP-seq_peaks.bed', dataset_dir(dataset))
write.table(as.data.frame(gr)[, 1:3], peak_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


# -----------------------------------------------------------------------------------
# [2019-10-07] 
# -----------------------------------------------------------------------------------
gs <- 'Etv2_MEF_reprogramming'; window_size <- 640; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'mm10'
source('analysis/seatac/windows.r'); gr <- load_windows(gs, window_size, bin_size, fs[1:2], fs[3])

G <- matrix(gr$num_reads, length(gr) / 5, 5)
i <- rowSums(G > 5) == 5
gr <- gr[rep(i, 5)]

Z <- gr$latent
k <- 4
u <- prcomp(Z)$x
set.seed(1); cls <- kmeans(Z, k)$cluster; mcols(gr)$cluster <- cls

H <- matrix(cls, length(gr) / 5, 5)

gr2 <- granges(gr[gr$group == 1])
mcols(gr2)$cluster <- H
gr_file <- sprintf('%s/Etv2_reprogramming.rds', dataset_dir('etv2_pioneer'))
saveRDS(gr2, gr_file)





