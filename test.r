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

latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 1; batch_effect <- FALSE
devtools::load_all('analysis/seatac/packages/seatac'); gr <- seatac(gr, latent_dim = latent_dim, n_components = n_components, prior = prior, epochs = epochs, batch_effect = batch_effect, beta = beta)
source('analysis/seatac/helper.r'); save_seatac_res(gr, gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)


# -----------------------------------------------------------------------------------
# [2019-07-19] t-SNE plot of the latent space of the windows
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox')
window_size <- 320; bin_size <- 10
latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
source('analysis/seatac/helper.r'); gr <- get_seatac_res(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)
Z <- mcols(gr)$latent[, , 1]
library(Rtsne); set.seed(1); y <- Rtsne(Z, check_duplicates = FALSE)$Y
Z <- mcols(gr)$latent[, , 1]
set.seed(1); nc <- 2; cls <- kmeans(Z, nc)$cluster
plot(y[, 1], y[, 2], bg = cls, pch = 21, col = cls, cex = 0.25, xaxt = 'n', yaxt = 'n')




# -----------------------------------------------------------------------------------
# [2019-06-26] Visualize the V-plot clusters for HMM prior
# -----------------------------------------------------------------------------------
library(gplots)
library(purrr)
#par(mfrow = c(10, 1), mar = c(0.25, 2, 0.25, 2))
par(mfrow = c(4, 6), mar = c(2, 4, 2, 4))
X <- mcols(gr)$counts %>% as.matrix()
dim(X) <- c(length(gr), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals, metadata(gr)$num_samples)
X <- aperm(X, c(1, 4, 2, 3))
dim(X) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window * metadata(gr)$n_intervals)
w <- 1 / rowSums(X); w[is.infinite(w)] <- 0
X <- as.matrix(Diagonal(x = w) %*% X)
dim(X) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window,  metadata(gr)$n_intervals)
Z <- mcols(gr)$latent
dim(Z) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$model$latent_dim * 2)
set.seed(1); nc <- 10; cls <- kmeans(Z, nc)$cluster
Y <- mcols(gr)$coverage
Y <- aperm(Y, c(1, 3, 2))
dim(Y) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window)
y <- c(50, 200, 400, 600, 670)
lapply(1:nc, function(k){
	Xm <- colSums(X[cls == k, , ], dims = 1)
	image(Xm, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = k)
	axis(2, at = (y - 50) / (670 - 50), label = y)
	axis(1, at = c(0, 0.5, 1))
	Ym <- colMeans(Y[cls == k, ], dims = 1)
	plot(Ym, type = 'b', main = k, lwd = 2)
})


table(mcols(gr)$cluster[, 1], mcols(gr)$cluster[, 2])
table(mcols(gr)$cluster[, 2], mcols(gr)$cluster[, 3])
table(mcols(gr)$cluster[, 3], mcols(gr)$cluster[, 4])



# -----------------------------------------------------------------------------------
# [2019-07-12] Compare the identified nucleosome and NFR regions with MNase-seq
# [2019-07-19] Including more epigenic data from MEF
# !!! Need to run on the lab queue
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox')
#gs <- c('Maza_MEF')
window_size <- 320; bin_size <- 10
latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
source('analysis/seatac/helper.r'); gr <- get_seatac_res(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)

devtools::load_all('packages/compbio'); bw_files <- get_bigwig_files()
library(ChIPpeakAnno); peaks <- reCenterPeaks(gr, width = 320)
extend <- 320; w <- 10
gss <- c('MEF_nucleosome', 'MEF_NoDox', 'MEF_MNase', 'MEF_Brg1', 'MEF_H3', 'MEF_H3K27me3', 'MEF_H3K27ac', 'MEF_H3K36me3', 'MEF_H3K9ac', 'MEF_H3K79me2', 'MEF_H3K4me2', 'MEF_H3K4me1', 'MEF_H3.3', 'MEF_P300', 'MEF_Dox_D1_Etv2')
#ps <- sprintf('%s_sample', metadata(gr)$peak_set); j <- sample(1:length(gr), 10000)
ps <- metadata(gr)$peak_set; j <- 1:length(gr)
devtools::load_all('packages/compbio'); n2m_files <- normalizeToMatrix_batch(peaks[j], peak_set = ps, bw_files[gss], extend = extend, w = w, mc.cores = 1, force = FALSE, target_ratio = 0.5)

library(EnrichedHeatmap)
library(circlize)
mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
cols <- c(
	'MEF_nucleosome' = 'blue',
	'MEF_NoDox' = 'blue',
	'MEF_MNase' = 'red',
	'MEF_H3' = 'red',
	'MEF_H3.3' = 'red',
	'MEF_Brg1' = 'red',
	'MEF_H3K9me3' = 'green',
	'MEF_H3K27me3' = 'green',
	'MEF_H3K27ac' = 'green',
	'MEF_H3K36me3' = 'green',
	'MEF_H3K9ac' = 'green',
	'MEF_H3K79me2' = 'green',
	'MEF_H3K4me2' = 'green',
	'MEF_H3K4me1' = 'green',
	'MEF_P300' = 'green',
	'MEF_Dox_D1_Etv2' = 'green'
)
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]], c(0, 0.99)), c('white', cols[g])))
names(col_fun) <- names(n2m_files)

Z <- mcols(peaks)$latent[j, , 1]
set.seed(1); nc <- 2; split_by <- kmeans(Z, nc)$cluster

#set.seed(1); i <- sample.int(nrow(Z), 20000)
set.seed(1); i <- 1:nrow(Z)
lgd <- Legend(at = c('C1', 'C2', 'C3'), title = "Clusters",  type = "lines", legend_gp = gpar(col = 2:4))
axis_name = c('-320', 'summit', '+320')
ta <- HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4), yaxis_side = 'right', yaxis_facing = 'inside'))
#gss <- c('MEF_MNase', 'MEF_NoDox', 'MEF_H3.3', 'MEF_H3', 'MEF_Brg1', 'MEF_Dox_D1_Etv2')
gss <- c('MEF_MNase', 'MEF_NoDox', 'MEF_Dox_D1_Etv2')
h <- Reduce('+', lapply(gss, function(s) EnrichedHeatmap(mat[[s]][i, ], col = col_fun[[s]], name = s, top_annotation = ta, axis_name = axis_name)))
draw(h, split = split_by[i], heatmap_legend_side = 'right', annotation_legend_list = list(lgd))



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

# -----------------------------------------------------------------------------------
# [2019-07-22] Use Seatac to two batch samples
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox', 'MEF_Dox_D1')
window_size <- 320; bin_size <- 10
source('analysis/seatac/helper.r'); gr <- get_fragment_size(gs, window_size = window_size, bin_size = bin_size, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, genome = BSgenome.Mmusculus.UCSC.mm10, exclude_exons = TRUE)

latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
devtools::load_all('analysis/seatac/packages/seatac'); gr <- seatac(gr[1:10000], latent_dim = latent_dim, n_components = n_components, prior = prior, epochs = epochs, batch_effect = batch_effect, beta = beta)
source('analysis/seatac/helper.r'); save_seatac_res(gr, gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)


window_size <- 320; bin_size <- 10
latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
source('analysis/seatac/helper.r'); gr <- get_seatac_res(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)


# -----------------------------------------------------------------------------------
# [2019-07-19] Use Seatac to analyze the Etv2 reprogramming
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox', 'MEF_Dox_D1', 'MEF_Dox_D2', 'MEF_Dox_D7', 'MEF_Dox_D7_Flk1pos')
window_size <- 320; bin_size <- 10
source('analysis/seatac/helper.r'); gr <- get_fragment_size(gs, window_size = window_size, bin_size = bin_size, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, genome = BSgenome.Mmusculus.UCSC.mm10, exclude_exons = TRUE)

latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
devtools::load_all('analysis/seatac/packages/seatac'); gr <- seatac(gr[1:10000], latent_dim = latent_dim, n_components = n_components, prior = prior, epochs = epochs, batch_effect = batch_effect, beta = beta)
source('analysis/seatac/helper.r'); save_seatac_res(gr, gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)


window_size <- 320; bin_size <- 10
latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
source('analysis/seatac/helper.r'); gr <- get_seatac_res(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)

