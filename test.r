library(roxygen2); library(devtools); devtools::create('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('packages/compbio')
    

# -----------------------------------------------------------------------------------
# [2019-06-24] VAE
# -----------------------------------------------------------------------------------

library(tensorflow)
tfe_enable_eager_execution(device_policy = 'silent')
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


# -----------------------------------------------------------------------------------
# [2019-07-11] segment the WT MEF
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox')
#gs <- c('Maza_MEF')
window_size <- 320; bin_size <- 10
source('analysis/seatac/helper.r'); gr <- get_fragment_size(gs, window_size = window_size, bin_size = bin_size, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, genome = BSgenome.Mmusculus.UCSC.mm10)

i <- 1000
plot(mcols(gr)$coverage[i, , 1], type = 'b'); image(matrix(mcols(gr)$counts[i, , ], 32, 32))

latent_dim <- 2; n_components <- 2; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
devtools::load_all('analysis/seatac/packages/seatac'); gr <- seatac(gr, latent_dim = latent_dim, n_components = n_components, prior = prior, epochs = epochs, batch_effect = batch_effect, beta = beta)
source('analysis/seatac/helper.r'); save_seatac_res(gr, gs, window_size, bin_size, min_reads_per_window, adjacent, latent_dim, n_components, prior, beta, epochs, batch_effect)



# -----------------------------------------------------------------------------------
# [2019-06-26] Visualize the V-plot clusters for HMM prior
# -----------------------------------------------------------------------------------
library(gplots)
library(purrr)
#par(mfrow = c(10, 1), mar = c(0.25, 2, 0.25, 2))
par(mfrow = c(2, 2), mar = c(2, 4, 2, 4))
X <- mcols(gr)$counts %>% as.matrix()
dim(X) <- c(length(gr), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals, metadata(gr)$num_samples)
X <- aperm(X, c(1, 4, 2, 3))
dim(X) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window * metadata(gr)$n_intervals)
w <- 1 / rowSums(X); w[is.infinite(w)] <- 0
X <- as.matrix(Diagonal(x = w) %*% X)
dim(X) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window,  metadata(gr)$n_intervals)
cls <- c(mcols(gr)$cluster)

y <- c(50, 200, 400, 600, 670)
lapply(1:metadata(gr)$model$n_components, function(k){
	Z <- colSums(X[cls == k, , ], dims = 1)
	image(Z, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = k)
	axis(2, at = (y - 50) / (670 - 50), label = y)
	axis(1, at = c(0, 0.5, 1))
})

table(mcols(gr)$cluster[, 1], mcols(gr)$cluster[, 2])
table(mcols(gr)$cluster[, 2], mcols(gr)$cluster[, 3])
table(mcols(gr)$cluster[, 3], mcols(gr)$cluster[, 4])



# -----------------------------------------------------------------------------------
# [2019-06-26] Visualize the V-plot clusters
# -----------------------------------------------------------------------------------
library(gplots)
#par(mfrow = c(10, 1), mar = c(0.25, 2, 0.25, 2))
par(mfrow = c(3, 5), mar = c(2, 0.25, 2, 0.25))
lapply(1:metadata(gr)$model$n_components, function(k) {
	X <- Reduce('+', lapply(1:metadata(gr)$num_samples, function(s){
		start <- 32 * 32 * (s - 1) + 1
		end <- 32 * 32 * s
		colSums(mcols(gr)$counts[mcols(gr)$cluster[, s] == k, start:end])
	}))
	X <- matrix(X, 32, 32)
	image(X, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
})
table(mcols(gr)$cluster)


# -----------------------------------------------------------------------------------
# [2019-07-12] Compare the identified nucleosome and NFR regions with MNase-seq
# Need to run on the lab queue
# -----------------------------------------------------------------------------------
#gs <- c('MEF_NoDox')
gs <- c('Maza_MEF')
window_size <- 320; bin_size <- 10; min_reads_per_window <- 20; adjacent <- 0
latent_dim <- 2; n_components <- 2; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
source('analysis/seatac/helper.r'); gr <- get_seatac_res(gs, window_size, bin_size, min_reads_per_window, adjacent, latent_dim, n_components, prior, beta, epochs, batch_effect)


source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
peaks <- reCenterPeaks(gr, width = 1)
#extend <- 320; w <- 10 
extend <- 1000; w <- 50 
gss <- c('MEF_NoDox', 'MEF_MNase')
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = metadata(gr)$peak_set, bw_files[gss], extend = extend, w = w, mc.cores = 1)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
cols <- c(
	'MEF_NoDox' = 'red',
	'MEF_MNase' = 'blue'
)
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]], c(0, 0.98)), c('white', cols[g])))
names(col_fun) <- names(n2m_files)

#i <- sample.int(nrow(mat[[1]]), 20000)
i <- 1:nrow(mat[[1]])
split_by <- mcols(peaks)$cluster[i, 1]

axis_name <- c('-1k', 'summit', '+1k')
#h <- EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], split = split_by, col = col_fun[['MEF_NoDox_ATAC']], name = 'MEF_NoDox_ATAC', axis_name = axis_name, pos_line = FALSE) +
#EnrichedHeatmap(mat[['MEF_MNase']][i, ], col = col_fun[['MEF_MNase']], name = 'MEF_MNase', axis_name = axis_name, pos_line = FALSE) 
h <- EnrichedHeatmap(mat[['MEF_MNase']][i, ], split = split_by, col = col_fun[['MEF_MNase']], name = 'MEF_MNase', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_NoDox']][i, ], col = col_fun[['MEF_NoDox']], name = 'MEF_NoDox_ATAC', axis_name = axis_name, pos_line = FALSE) 
draw(h, heatmap_legend_side = 'right')


M <- as.matrix(mat[['MEF_MNase']])
plot(colMeans(M[mcols(peaks)$cluster[, 1] == 1, ]), type = 'l', lwd = 2); abline(v = 20.5, lwd = 2, lty = 2)

plot(colMeans(M[mcols(peaks)$cluster[, 1] == 2, ]), type = 'l', lwd = 2); abline(v = 20.5, lwd = 2, lty = 2)




# -----------------------------------------------------------------------------------
# [2019-06-26] Prepare the public available ATAC-seq data for MEF and mESC
# -----------------------------------------------------------------------------------
dataset <- 'dataset=Maza_version=20170302a'
source('sra.r'); d <- read.dataset(dataset, touch = TRUE)
d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
treatment_files <- d[d$group %in% c('mMEF'), 'bam.file']

merged_treatment_file <- 'analysis/seatac/data/Maza_MEF.bam'
source('chipseq.r'); merge.bam.files(treatment_files, merged_treatment_file)	# index will be created automatically

source('analysis/seatac/helper.r'); base.name <- 'analysis/seatac/data/Maza_MEF'
source('chipseq.r');  macs2.callpeak(merged_treatment_file, base.name, format = 'BAMPE', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = TRUE,  shift = -100, extsize = 200)



