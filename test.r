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
# [2019-07-11] Training VAE on MEF ATAC-seq data
# -----------------------------------------------------------------------------------
#gs <- 'MEF_NoDox'; ps <- 'MEF_NoDox'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 2000; window_size <- 400; step_size <- 200
#gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 400; window_size <- 400; bin_size <- 10; step_size <- 400; latent_dim <- 2; n_components <- 20; epochs <- 100; batch_effect <- FALSE
#gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 400; window_size <- 400; bin_size <- 5; step_size <- 400; latent_dim <- 2; n_components <- 20; epochs <- 100; batch_effect <- FALSE
gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS_0_400'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 400; window_size <- 400; bin_size <- 10; step_size <- NA; latent_dim <- 2; n_components <- 20; epochs <- 100; batch_effect <- FALSE
#gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 1600; window_size <- 400; step_size <- 400
#gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 800; window_size <- 800; step_size <- 400; latent_dim <- 2; n_components <- 20; epochs <- 100; batch_effect <- FALSE
#gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS'; min_reads_per_window <- 15; min_reads_coverage <- 30; expand <- 800; window_size <- 800; step_size <- 400; latent_dim <- 2; n_components <- 20; epochs <- 100; batch_effect <- FALSE
#gs <- 'MEF_NoDox'; ps <- 'MEF_active_TSS2'; min_reads_per_window <- 50; min_reads_coverage <- 30; expand <- 2000

source('analysis/seatac/helper.r'); peaks <- read_peaks(ps)
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, peaks, genome = BSgenome.Mmusculus.UCSC.mm10, expand = expand * 2)
source('analysis/seatac/helper.r'); windows <- read_windows(gs, ga, peaks, expand = expand, bin_size = bin_size, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, min_reads_per_window = min_reads_per_window)
source('analysis/seatac/helper.r'); model_dir <- model_dir_name(gs, ps, expand, latent_dim, n_components, batch_effect, window_size, step_size, min_reads_per_window, min_reads_coverage)

# run the model 
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, n_components = n_components, window_size = window_size, step_size = step_size, epochs = epochs, batch_effect = batch_effect, min_reads_per_window = min_reads_per_window, min_reads_coverage = min_reads_coverage)
devtools::load_all('analysis/seatac/packages/seatac'); saveModel(model, model_dir)

# load the model
devtools::load_all('analysis/seatac/packages/seatac'); model <- loadModel(model_dir)


# -----------------------------------------------------------------------------------
# [2019-07-26] predict and segment by HMM
# -----------------------------------------------------------------------------------
step_size <- 50
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows[1:200], window_size = window_size, step_size = step_size)

devtools::load_all('analysis/seatac/packages/seatac'); gr <- gr %>% segment(k = 2)



# -----------------------------------------------------------------------------------
# [2019-07-25] Predicting on the entire MEF peaks region
# -----------------------------------------------------------------------------------
step_size <- 50
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows, window_size = window_size, step_size = step_size)
gr_file <- sprintf('%s/gr_step_size=%d.rds', model_dir, step_size)
flog.info(sprintf('writing %s', gr_file))
saveRDS(gr, gr_file)


par(mfcol = c(4, 2))
i <- 1
yy <- c(50, 200, 400, 600, 670)
image(matrix(as.matrix(mcols(gr)$counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
image(matrix(as.matrix(mcols(gr)$fitted_counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
plot(mcols(gr)$coverage[i,], type = 'b', lwd = 2)
plot(mcols(gr)$fitted_coverage[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")



# -----------------------------------------------------------------------------------
# [2019-07-25] Predict and visualize V plot on a single genomic region
# -----------------------------------------------------------------------------------
source('analysis/seatac/helper.r')
gr <- GRanges(seqnames = 'chr3', range = IRanges(start = 137269432, end = 137297091))
expand <- 10000
gr <- resize(gr, fix = 'center', width = expand)
devtools::load_all('analysis/seatac/packages/seatac'); ga <- readBAM(list_files()['MEF_NoDox', 'bam'], gr, genome = BSgenome.Mmusculus.UCSC.mm10, expand = expand + 5000)
devtools::load_all('analysis/seatac/packages/seatac');  gr <- readFragmentSizeMatrix(ga, gr, window_size = expand)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(gr, window_size = 400, step_size = 50)

devtools::load_all('packages/compbio'); bw_files <- get_bigwig_files()
bw_file <- bw_files['MEF_MNase']
cvg <- rtracklayer::import(bw_file, which = gr, as = 'Rle')
seqinfo(gr) <- seqinfo(cvg)
x_mnase <- mean(cvg[unlist(slidingWindows(gr, width = 10, step = 10))])


par(mfcol = c(5, 1))
i <- 1
yy <- c(50, 200, 400, 600, 670)
image(matrix(as.matrix(mcols(gr)$counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
image(matrix(as.matrix(mcols(gr)$fitted_counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
plot(mcols(gr)$coverage[i,], type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(mcols(gr)$fitted_coverage[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(x_mnase, type = 'b', lwd = 2, xaxs="i", yaxs="i")



# -----------------------------------------------------------------------------------
# [2019-07-25] Compare MEF and D1 Etv2 V-plot for a specific genomic region
# -----------------------------------------------------------------------------------
source('analysis/seatac/helper.r')
gr <- GRanges(seqnames = 'chr3', range = IRanges(start = 137269432, end = 137297091))
expand <- 10000
gr <- resize(gr, fix = 'center', width = expand)
devtools::load_all('analysis/seatac/packages/seatac'); ga <- readBAM(list_files()['MEF_NoDox', 'bam'], gr, genome = BSgenome.Mmusculus.UCSC.mm10, expand = expand + 5000)
devtools::load_all('analysis/seatac/packages/seatac');  gr <- readFragmentSizeMatrix(ga, gr, window_size = expand)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(gr, window_size = 400, step_size = 50)
devtools::load_all('analysis/seatac/packages/seatac'); ga2 <- readBAM(list_files()['MEF_Dox_D1', 'bam'], gr, genome = BSgenome.Mmusculus.UCSC.mm10, expand = expand + 5000)
devtools::load_all('analysis/seatac/packages/seatac');  gr2 <- readFragmentSizeMatrix(ga2, gr, window_size = expand)
devtools::load_all('analysis/seatac/packages/seatac'); gr2 <- model %>% predict(gr2, window_size = 400, step_size = 50)

par(mfcol = c(4, 1))
i <- 1
yy <- c(50, 200, 400, 600, 670)
#image(matrix(as.matrix(mcols(gr)$counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
#axis(2, at = (yy - 50) / (670 - 50), label = yy)
#axis(1, at = c(0, 0.5, 1))
image(matrix(as.matrix(mcols(gr)$fitted_counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))

image(matrix(as.matrix(mcols(gr2)$fitted_counts[i, ]), metadata(gr2)$n_bins_per_window, metadata(gr2)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))

#plot(mcols(gr)$coverage[i,], type = 'b', lwd = 2)
plot(mcols(gr)$fitted_coverage[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")

#image(matrix(as.matrix(mcols(gr2)$counts[i, ]), metadata(gr2)$n_bins_per_window, metadata(gr2)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
#axis(2, at = (yy - 50) / (670 - 50), label = yy)
#axis(1, at = c(0, 0.5, 1))
#plot(mcols(gr2)$coverage[i,], type = 'b', lwd = 2)
plot(mcols(gr2)$fitted_coverage[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")



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


															#     y <- model$decoder$coverage(Zb[i, , drop = FALSE])$sample(1000L) %>% as.array()
															#     plot(drop(colSums(y, dim = 1)), type = 'b', lwd = 2)


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
bin_size <- 10; expand <- 1000
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, genome = BSgenome.Mmusculus.UCSC.mm10)
source('analysis/seatac/helper.r'); windows <- read_windows(gs, ga, window_size = expand, bin_size = bin_size, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, exclude_exons = TRUE)

latent_dim <- 2; n_components <- 20; epochs <- 100; batch_effect <- FALSE; window_size <- 400; step_size <- 200
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows[1:10000], latent_dim = latent_dim, n_components = n_components, window_size = window_size, step_size = step_size, epochs = epochs, batch_effect = batch_effect)

devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows, step_size = step_size)


source('analysis/seatac/helper.r'); save_seatac_res(gr, gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)


window_size <- 320; bin_size <- 10
latent_dim <- 2; n_components <- 10; prior <- 'gmm'; beta <- 1; epochs <- 100; batch_effect <- FALSE
source('analysis/seatac/helper.r'); gr <- get_seatac_res(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect)


# -----------------------------------------------------------------------------------
# [2019-08-01] Look at the V-plot at the TSS region
# -----------------------------------------------------------------------------------
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows, window_size = window_size, step_size = step_size)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- gr %>% segment(k = 2)

minus <- as.logical(strand(gr) == '-')
X <- array(as.matrix(mcols(gr)$counts), dim = c(length(gr), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals))
X[minus, , ] <- X[minus, metadata(gr)$n_bins_per_window:1, ]
Xp <- array(as.matrix(mcols(gr)$fitted_counts), dim = c(length(gr), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals))
Xp[minus, , ] <- Xp[minus, metadata(gr)$n_bins_per_window:1, ]


image(X[i, , ], col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
abline(v = 0.5)
image(Xp[i, , ], col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
abline(v = 0.5)
gr[i]

plot(S[i, ])
abline(v = 0.5)



j <- mcols(gr)$expression[, 1] > 1
image(colMeans(Xp[j, , ], 1), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
abline(v = 0.5)
image(colMeans(X[j, , ], 1), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
abline(v = 0.5)



# -----------------------------------------------------------------------------------
# [2019-08-01] Look at the dinucleosome distribution of identified nucleosome
# -----------------------------------------------------------------------------------
step_size <- 50
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows[1:1000], window_size = window_size, step_size = step_size)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- gr %>% segment(k = 2)
devtools::load_all('analysis/seatac/packages/seatac'); gr2 <- gr %>% call_nucleosome()



