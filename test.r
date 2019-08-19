library(roxygen2); library(devtools); devtools::create('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('packages/compbio')
    

# -----------------------------------------------------------------------------------
# [2019-06-24] VAE
# -----------------------------------------------------------------------------------

library(tensorflow)
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)


# -----------------------------------------------------------------------------------
# [2019-07-11] Training nucleosome classifier on mESC ATAC-seq data and chemicall mapping
# data
# -----------------------------------------------------------------------------------
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC_chr1-3'; window_size <- 320; bin_size <- 10; step_size <- 20; mr <- 5; mc <- 0; bs <- 256; ns <- 1; seed <- 1
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC_chr1'; window_size <- 320; bin_size <- 5; step_size <- 40; mr <- 5; mc <- 0; bs <- 128
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC_chr1-3'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1; cw <- 20
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1; cw <- 10
#gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1; cw <- 147
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1; cw <- 5

latent_dim <- 2; epochs <- 30

source('analysis/seatac/helper.r'); peaks <- read_peaks(ps)
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
source('analysis/seatac/helper.r'); windows <- prepare_training_windows(peaks, expand = expand, window_size = window_size, core_width = cw, negative_sample_ratio = ns)
source('analysis/seatac/helper.r'); windows <- readFragmentSizeMatrix(ga, windows, window_size = window_size, bin_size = bin_size)
windows <- windows[mcols(windows)$num_reads >= mr & mcols(windows)$mean_coverage >= mc]
source('analysis/seatac/helper.r'); model_dir <- model_dir_name(gs, ps, latent_dim, expand, window_size, mr, mc, bin_size, ns, cw)

# run the model 
devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, epochs = epochs, batch_size = bs)
devtools::load_all('analysis/seatac/packages/seatac'); save_model(model, model_dir)


# load the model
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)


devtools::load_all('analysis/seatac/packages/seatac'); model %>% predict(windows[1:100])

# -----------------------------------------------------------------------------------
# [2019-08-16] Predicting the nucleosome positions ATAC-seq peak region on mESC
# Compare the ROC curve  between seatac and NucleoATAC
# -----------------------------------------------------------------------------------

# --- prepare all the positive and negative samples
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; bs <- 256; ns <- 1
latent_dim <- 2; epochs <- 20
source('analysis/seatac/helper.r'); peaks <- read_peaks(ps)
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
source('analysis/seatac/helper.r'); windows <- prepare_training_windows(peaks, expand = expand, window_size = window_size, negative_sample_ratio = ns)
source('analysis/seatac/helper.r'); windows <- readFragmentSizeMatrix(ga, windows, window_size = window_size, bin_size = bin_size)

# --- train the model based on each of the 4/5 the data and predict the nucleosome score on the remaining 1/5
cross <- 5
set.seed(1); groups <- sample(1:cross, length(windows), replace = TRUE)
wins <- lapply(1:cross, function(i){
	windows_train <- windows[groups != i]
	windows_test <- windows[groups == i]
	windows_train <- windows_train[mcols(windows_train)$num_reads >= mr & mcols(windows_train)$mean_coverage >= mc]
	model <- seatac(windows_train, latent_dim = latent_dim, epochs = epochs, batch_size = bs)
	windows_test <- model %>% predict(windows_test)
})
wins <- Reduce('c', wins)
wins <- add.seqinfo(wins, 'mm10')
devtools::load_all('packages/compbio'); win_file <- sprintf('%s/Maza_mESC_seatac_5fold_cv.rds', project_dir('seatac'))
saveRDS(wins, win_file)

# --- load the predicted nucleosome scores
devtools::load_all('packages/compbio'); win_file <- sprintf('%s/Maza_mESC_seatac_5fold_cv.rds', project_dir('seatac'))
wins <- readRDS(win_file)

# --- read ~80M called nucleosomes by the chemical mapping method
nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))  # a simple seqnames/start/end format
nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
nuc <- add.seqinfo(nuc, 'mm10')

# --- read the NCP scores over the windows
ncp_file <- sprintf('%s/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz', sra.run.dir('GSM2183909'))
flog.info(sprintf('reading %s', ncp_file))
ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(resize(wins, width = window_size, fix = 'center')))
ncp <- resize(ncp, width = 1, fix = 'center')
ncp <- add.seqinfo(ncp, 'mm10')
cvg_ncp <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))

# --- read the nucleosome signal by NucleoATAC
dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
na <- add.seqinfo(na, 'mm10')
cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))


# --- Compar the sensitivity / specificity 
wins <- resize(wins, fix = 'center', width = 147)
mcols(wins)$label <- wins %over% nuc	# the ground truth
mcols(wins)$label_nucleoatac <- wins %over% na	# whether marked as nucleosome by NucleoATAC
mcols(wins)$nucleoatac_score <- mean(cvg_nucleoatac[wins])	# the continuous call by NucleoATAC
mcols(wins)$ncp_score <- mean(cvg_ncp[wins])

res <- do.call('rbind', lapply(seq(0.01, 0.99, by = 0.01), function(h){
	y_true <- as.numeric(mcols(windows_test)$label)
	Y <- data.frame(
		nucleoatac = mcols(windows_test)$nucleoatac_score > quantile(mcols(windows_test)$nucleoatac_score, h),
		seatac = mcols(windows_test)$nucleosome_score> quantile(mcols(windows_test)$nucleosome_score, h)
	)
	cbind(h = h, do.call('rbind', lapply(1:ncol(Y), function(i){
		y_pred <- Y[, i]
		tp <- sum(y_true & y_pred)
		fp <- sum(!y_true & y_pred)
		fn <- sum(y_true & !y_pred)
		tn <- sum(!y_true & !y_pred)
		sensitivity <- tp / (tp + fn)
		specificity <- tn / (tn + fp)
		data.frame(n = sum(y_pred), pct_true = mean(y_pred), method = colnames(Y)[i], sensitivity = sensitivity, specificity = specificity, ppv = tp / (tp + fp), f1 = 2 * tp / (2 * tp + fp + fn))
	})))
}))

plot(1 - res$specificity[res$method == 'seatac'], res$sensitivity[res$method == 'seatac'], type = 'l', lwd = 2, col = 'blue', xlim = c(0, 1), ylim = c(0, 1))
lines(1 - res$specificity[res$method == 'nucleoatac'], res$sensitivity[res$method == 'nucleoatac'], lwd = 2, col = 'black')

plot(res$pct_true[res$method == 'seatac'], res$ppv[res$method == 'seatac'], type = 'l', lwd = 2, col = 'blue', xlim = c(0, 1), ylim = c(0.2, 1))
lines(res$pct_true[res$method == 'nucleoatac'], res$ppv[res$method == 'nucleoatac'], type = 'l', lwd = 2, col = 'black')


# -----------------------------------------------------------------------------------
# [2019-08-18] Predict the nucleosome score on every 10 bp bins of all
# 320-bp windows surrounding each ATAC-seq summits
# -----------------------------------------------------------------------------------
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand0 <- 320; window_size <- 320; bin_size <- 10; step_size <- 10
source('analysis/seatac/helper.r'); peaks <- read_peaks(ps)
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, genome = BSgenome.Mmusculus.UCSC.mm10)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- read_windows(ga, peaks, expand = expand0, window_size = window_size, step_size = step_size, bin_size = bin_size)

gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; latent_dim <- 2; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; ns <- 1; cw <- 5
source('analysis/seatac/helper.r'); model_dir <- model_dir_name(gs, ps, latent_dim, expand, window_size, mr, mc, bin_size, ns, cw)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
windows <- model %>% predict(windows)
windows <- add.seqinfo(windows, 'mm10')

#--- Add the mean NCP score for each bin
ncp_file <- sprintf('%s/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz', sra.run.dir('GSM2183909'))
flog.info(sprintf('reading %s', ncp_file))
ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(resize(peaks, width = window_size, fix = 'center')))
ncp <- resize(ncp, width = 1, fix = 'center')
ncp <- add.seqinfo(ncp, 'mm10')
cvg_ncp <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
mcols(windows)$ncp_score <- mean(cvg_ncp[resize(windows, fix = 'center', width = step_size)])

# --- Add the known nucleosome indicator (ground truth)
nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))  # a simple seqnames/start/end format
nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
nuc <- add.seqinfo(nuc, 'mm10')
mcols(windows)$label <- resize(windows, fix = 'center', width = step_size) %over% nuc

# --- Add nucleosome signal by NucleoATAC
dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
na <- add.seqinfo(na, 'mm10')
cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))
mcols(windows)$score_nucleoatac <- mean(cvg_nucleoatac[resize(windows, fix = 'center', width = step_size)])

# --- Add the nucleosome identified by NuceloATAC
dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.occpeaks.bed.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
occ <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
occ <- GRanges(seqnames = occ[, 1], range = IRanges(occ[, 2], occ[, 3]), score = occ[, 4])
occ <- occ[occ %over% resize(peaks, width = window_size, fix = 'center')]
occ <- add.seqinfo(occ, 'mm10')
mcols(windows)$label_nucleoatac <- resize(windows, fix = 'center', width = step_size) %over% occ


# --- look at the mean NCP score of positive/negative predicted bins
label_pred2 <- mcols(windows)$nucleosome_score > quantile(mcols(windows)$nucleosome_score, 0.95)
boxplot(list(mcols(windows)$ncp_score[label_pred2], mcols(windows)$ncp_score[!label_pred2]), ylim = c(0, 0.25), outline = FALSE)

boxplot(list(mcols(windows)$ncp_score[mcols(windows)$label_nucleoatac], mcols(windows)$ncp_score[!mcols(windows)$label_nucleoatac]), ylim = c(0, 0.25), outline = FALSE)
boxplot(list(mcols(windows)$ncp_score[mcols(windows)$label], mcols(windows)$ncp_score[!mcols(windows)$label]), ylim = c(0, 0.25), outline = FALSE)

# --- ROC and PPV 
j <- mcols(windows)$mean_coverage > 5 & mcols(windows)$num_reads > 5
res <- do.call('rbind', lapply(seq(0.9, 0.999, by = 0.001), function(h){
	y_true <- as.numeric(mcols(windows)$label)
	Y <- data.frame(
		nucleoatac = mcols(windows)$score_nucleoatac > quantile(mcols(windows)$score_nucleoatac, h),
		seatac = mcols(windows)$nucleosome_score> quantile(mcols(windows)$nucleosome_score, h)
	)
	y_true <- y_true[j]; Y <- Y[j, ]
	cbind(h = h, do.call('rbind', lapply(1:ncol(Y), function(i){
		y_pred <- Y[, i]
		tp <- sum(y_true & y_pred)
		fp <- sum(!y_true & y_pred)
		fn <- sum(y_true & !y_pred)
		tn <- sum(!y_true & !y_pred)
		sensitivity <- tp / (tp + fn)
		specificity <- tn / (tn + fp)
		data.frame(n = sum(y_pred), pct_true = mean(y_pred), method = colnames(Y)[i], sensitivity = sensitivity, specificity = specificity, ppv = tp / (tp + fp), f1 = 2 * tp / (2 * tp + fp + fn))
	})))
}))


plot(1 - res$specificity[res$method == 'seatac'], res$sensitivity[res$method == 'seatac'], type = 'l', lwd = 2, col = 'blue', xlim = c(0, 1), ylim = c(0, 1))
lines(1 - res$specificity[res$method == 'nucleoatac'], res$sensitivity[res$method == 'nucleoatac'], lwd = 2, col = 'black')

plot(res$pct_true[res$method == 'seatac'], res$ppv[res$method == 'seatac'], type = 'l', lwd = 2, col = 'blue', xlim = c(0, 0.1))
lines(res$pct_true[res$method == 'nucleoatac'], res$ppv[res$method == 'nucleoatac'], type = 'l', lwd = 2, col = 'black')



# -----------------------------------------------------------------------------------
# [2019-08-18] Look at the AT/TT/TA pattern surrounding the called nucleosome center
# -----------------------------------------------------------------------------------
gr <- list()

label_pred2 <- mcols(windows)$nucleosome_score > quantile(mcols(windows)$nucleosome_score, 0.90)
gr[['seatac']] <- windows[label_pred2]

dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.occpeaks.bed.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
occ <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
occ <- GRanges(seqnames = occ[, 1], range = IRanges(occ[, 2], occ[, 3]), score = occ[, 4])
gr[['nucleoatac']] <- occ[occ %over% resize(peaks, width = window_size, fix = 'center')]
gr[['nucleoatac']] <- add.seqinfo(gr[['nucleoatac']], 'mm10')

nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))  # a simple seqnames/start/end format
nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
gr[['chemical']] <- nuc[nuc %over% resize(peaks, width = window_size, fix = 'center')]
gr[['chemical']] <- add.seqinfo(gr[['chemical']], 'mm10')




# -----------------------------------------------------------------------------------
# [2019-08-18] read the NCP scores over the *peaks*
# -----------------------------------------------------------------------------------

w <- 147
dinuc <- expand.grid(c('A', 'C', 'G', 'T'), c('A', 'C', 'G', 'T'))
dinuc <- sprintf('%s%s', dinuc[, 1], dinuc[, 2])
Y <- list()

for (g in c('nucleoatac', 'chemical', 'seatac')){
	mcols(gr[[g]])$sequence <- getSeq(Mmusculus, resize(gr[[g]], fix = 'center', width = w))
	Y[[g]] <- do.call('cbind', lapply(1:(w - 1), function(i){
		table(factor(as.character(subseq(mcols(gr[[g]])$sequence, start = i, width = 2)), dinuc)) / length(gr[[g]])
	}))
}

WW <- do.call('rbind', lapply(1:length(Y), function(i){
	ww <- colSums(Y[[i]][c('AA', 'AT', 'TA', 'TT'), ])
#	ww <- (ww - min(ww)) / (max(ww) - min(ww))
	ww
}))

plot(WW[1, ], type = 'b')
plot(WW[2, ], type = 'b')
plot(WW[3, ], type = 'b')


# -----------------------------------------------------------------------------------
# [2019-08-18] Check see if there is bias based on AT at called nucleosome
# -----------------------------------------------------------------------------------
library(seqinr)
w2 <- 10
mcols(windows)$sequence <- getSeq(Mmusculus, resize(windows, fix = 'center', width = w2))
mcols(windows)$gc <- rowSums(alphabetFrequency(mcols(windows)$sequence)[, 2:3]) / w2

boxplot(split(mcols(windows)$gc, list(mcols(windows)$label)))
boxplot(split(mcols(windows)$gc, list(label_pred2)))
boxplot(split(mcols(windows)$gc, list(mcols(windows)$label_nucleoatac)))





# -----------------------------------------------------------------------------------
# [2019-08-18] Look at the called nucleosome positions surrounding the TSS
# -----------------------------------------------------------------------------------
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; latent_dim <- 2; expand0 <- 500; window_size <- 320; bin_size <- 10; bs <- 256; step_size <- 5
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
pmt <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream = expand0 / 2, downstream = expand0 / 2)
source('analysis/seatac/helper.r'); peaks <- pmt[pmt %over% read_peaks(ps)]
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, genome = BSgenome.Mmusculus.UCSC.mm10)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- read_windows(ga, peaks, expand = expand0, window_size = window_size, step_size = step_size, bin_size = bin_size)

# --- load the trained model on mESC and predict the nucleosome score
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; latent_dim <- 2; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; ns <- 1; cw <- 147
source('analysis/seatac/helper.r'); model_dir <- model_dir_name(gs, ps, latent_dim, expand, window_size, mr, mc, bin_size, ns, cw)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)
windows <- add.seqinfo(windows, 'mm10')

X <- matrix(mcols(windows)$nucleosome_score, length(peaks), expand0 / step_size, byrow = TRUE)
X <- matrix(rowMeans(mcols(windows)$vplot_center[, 1:7]), length(peaks), expand0 / step_size, byrow = TRUE)

# --- read the NCP score
ncp_file <- sprintf('%s/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz', sra.run.dir('GSM2183909'))
flog.info(sprintf('reading %s', ncp_file))
ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(resize(peaks, width = expand0 + 20, fix = 'center')))
ncp <- resize(ncp, width = 1, fix = 'center')
ncp <- add.seqinfo(ncp, 'mm10')
cvg_ncp <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
X1 <- as(as(cvg_ncp[resize(peaks, fix = 'center', width = expand0)], 'RleViews'), 'matrix')

# --- read the nucleosome signal by NucleoATAC
dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
na <- add.seqinfo(na, 'mm10')
cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))
X2 <- as(as(cvg_nucleoatac[resize(peaks, fix = 'center', width = expand0)], 'RleViews'), 'matrix')

plot(colMeans(X[as.character(strand(peaks)) == '-', ]), type = 'b'); abline(v = expand0 / step_size / 2)
plot(colMeans(X1[as.character(strand(peaks)) == '-', ]), type = 'b'); abline(v = expand0 / 2)
plot(colMeans(X2[as.character(strand(peaks)) == '-', ]), type = 'b'); abline(v = expand0 / 2)



# -----------------------------------------------------------------------------------
# [2019-08-18] Look at continuous score surrounding each known nucleosome
# -----------------------------------------------------------------------------------
nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))  # a simple seqnames/start/end format
nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
nuc <- add.seqinfo(nuc, 'mm10')
ps <- 'Maza_mESC'; window_size <- 320
source('analysis/seatac/helper.r'); peaks <- nuc[nuc %within% resize(read_peaks(ps), fix = 'center', width = window_size)]

# --- score by seatac
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; expand0 <- 500; window_size <- 320; bin_size <- 10; step_size <- 20
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, genome = BSgenome.Mmusculus.UCSC.mm10)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- read_windows(ga, peaks, expand = expand0, window_size = window_size, step_size = step_size, bin_size = bin_size)
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; latent_dim <- 2; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; ns <- 1; cw <- 147
source('analysis/seatac/helper.r'); model_dir <- model_dir_name(gs, ps, latent_dim, expand, window_size, mr, mc, bin_size, ns, cw)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- model %>% predict(windows)
windows <- add.seqinfo(windows, 'mm10')
X <- matrix(mcols(windows)$nucleosome_score > 0, length(peaks), expand0 / step_size, byrow = TRUE)


# --- read the nucleosome signal by NucleoATAC
dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
na <- add.seqinfo(na, 'mm10')
cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))
X2 <- as(as(cvg_nucleoatac[resize(peaks, fix = 'center', width = expand0)], 'RleViews'), 'matrix')


mcols(windows)$label <- resize(windows, fix = 'center', width = 147) %over% nuc
mcols(windows)$label_nucleoatac <- resize(windows, fix = 'center', width = 20) %over% occ

mcols(windows)$ncp_score <- mean(cvg_ncp[resize(windows, fix = 'center', width = 20)])
mcols(windows)$nucleoatac_score <- mean(cvg_nucleoatac[resize(windows, fix = 'center', width = 147)])	# the continuous call by NucleoATAC






# -----------------------------------------------------------------------------------
# [2019-07-25] Predicting on the all MEF peaks 
# -----------------------------------------------------------------------------------
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% encode(windows, window_size = window_size, step_size = step_size)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% decode(gr)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- gr %>% segment()


center <- 31:33
is_nucleosome <- rowSums(mcols(gr)$state[, center, drop = FALSE] == 1) == 3
is_nfr <- rowSums(mcols(gr)$state[, center, drop = FALSE] == -1) == 3
image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nfr, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)


image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nucleosome, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nfr, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)

w <- width(gr)[1]
s <- getSeq(Mmusculus, gr)
dinuc <- expand.grid(c('A', 'C', 'G', 'T'), c('A', 'C', 'G', 'T'))
dinuc <- sprintf('%s%s', dinuc[, 1], dinuc[, 2])

Y <- do.call('cbind', lapply(1:(w - 1), function(i){
	table(factor(as.character(subseq(s[is_nucleosome, ], start = i, width = 2)), dinuc)) / sum(is_nucleosome)
}))
ww <- colSums(Y[c('AA', 'AT', 'TA', 'TT'), ])
ww <- (ww - min(ww)) / (max(ww) - min(ww))
rr <- colSums(Y[c('GG', 'GC', 'CG', 'CC'), ])
rr <- (rr - min(rr)) / (max(rr) - min(rr))

plot(ww, type = 'l', col = 'red')
lines(rr, type = 'l', col = 'blue')
abline(v = 160, lwd = 2, lty = 2

plot(colSums(Y[c('AA', 'AT', 'TA', 'TT'), ]), type = 'b', col = 'red')


image(matrix(colMeans(as.matrix(mcols(gr)$counts[clust == h, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = h)



# -----------------------------------------------------------------------------------
# [2019-07-12] Compare the identified nucleosome and NFR regions with MNase-seq
# [2019-07-19] Including more epigenic data from MEF
# [2019-08-05] 
# -----------------------------------------------------------------------------------
devtools::load_all('packages/compbio'); bw_files <- get_bigwig_files()
library(ChIPpeakAnno); peaks <- reCenterPeaks(gr, width = 320)
mcols(peaks)$is_nucleosome <- is_nucleosome
mcols(peaks)$is_nfr <- is_nfr
peaks <- peaks[mcols(peaks)$is_nucleosome | mcols(peaks)$is_nfr]

extend <- 320; w <- 10
gss <- c('MEF_nucleosome', 'MEF_NoDox', 'MEF_MNase', 'MEF_H3', 'MEF_H3.3')
devtools::load_all('packages/compbio'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'MEF_nucloesome_or_nfr', bw_files[gss], extend = extend, w = w, mc.cores = 1, force = FALSE, target_ratio = 0.5)

library(EnrichedHeatmap)
library(circlize)
mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
cols <- c(
	'MEF_MNase' = 'red',
	'MEF_H3' = 'red',
	'MEF_NoDox' = 'blue'
)
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]], c(0, 0.99)), c('white', cols[g])))
names(col_fun) <- names(n2m_files)


lgd <- Legend(at = c('FALSE', 'TRUE'), title = 'Clusters',  type = 'lines', legend_gp = gpar(col = 2:4))
axis_name = c('-320', 'summit', '+320')
ta <- HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c('black', 'red')), yaxis_side = 'right', yaxis_facing = 'inside'))
gss <- c('MEF_MNase', 'MEF_H3', 'MEF_NoDox')
h <- Reduce('+', lapply(gss, function(s) EnrichedHeatmap(mat[[s]], col = col_fun[[s]], name = s, top_annotation = ta, axis_name = axis_name)))
draw(h, split = factor(mcols(peaks)$is_nucleosome), heatmap_legend_side = 'right', annotation_legend_list = list(lgd))

barplot(c(7910, 143))


par(mfcol = c(6, 6))
i <- 2
yy <- c(50, 200, 400, 600, 670)
image(matrix(as.matrix(mcols(gr)$counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
image(matrix(as.matrix(mcols(gr)$fitted_counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
plot(mcols(gr)$coverage[i,], type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(mcols(gr)$fitted_coverage[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(mcols(gr)$state[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i", ylim = c(-1, 1))
plot(mcols(gr)$z_score[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")




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


par(mfcol = c(4, 6))
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


# -----------------------------------------------------------------------------------
# [2019-07-25] Look at the changes between NFR and nucleosome during Etv2 reprogramming
# -----------------------------------------------------------------------------------
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% encode(windows, window_size = window_size, step_size = step_size, batch_size = 2^11)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% decode(gr, batch_size = 2^11)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- gr %>% segment()

centers <- (32 - 1):(32 + 1)
is_nucleosome <- rowSums(mcols(gr)$z_score[, centers, drop = FALSE] > 1) == length(centers)
is_nfr <- rowSums(mcols(gr)$z_score[, centers, drop = FALSE] < -1) == length(centers)
#is_nucleosome <- rowSums(mcols(gr)$state[, centers, drop = FALSE] == 1) == length(centers)
#is_nfr <- rowSums(mcols(gr)$state[, centers, drop = FALSE] == -1) == length(centers)
table(is_nucleosome)
table(is_nfr)


image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nucleosome, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nfr, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)


s <- rep('unknown', length(gr))
s[is_nucleosome] <- 'nucleosome'
s[is_nfr] <- 'nfr'
S <- matrix(s, length(gr) / metadata(gr)$num_samples, metadata(gr)$num_samples)
S <- S[, c(1, 2)]
S <- S[rowSums(S == 'unknown') == 0, ]
#y <- sprintf('%s-%s-%s', S[, 1], S[, 2], S[, 3])
y <- sprintf('%s-%s', S[, 1], S[, 2])
table(y)



table(s[mcols(gr)$group == 1], s[mcols(gr)$group == 2])
table(s[mcols(gr)$group == 2], s[mcols(gr)$group == 3])
table(s[mcols(gr)$group == 3], s[mcols(gr)$group == 4])
table(s[mcols(gr)$group == 1], s[mcols(gr)$group == 5])

image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nfr, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)




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
# [2019-08-02] Prepare the 
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox', 'MEF_Dox_D1', 'MEF_Dox_D2', 'MEF_Dox_D7', 'MEF_Dox_D7_Flk1pos'); ps <- c('MEF_Dox_D1')
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
# [2019-08-02] Look at clusters of V-plot
# -----------------------------------------------------------------------------------
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% encode(windows, window_size = window_size, step_size = step_size)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% decode(gr)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- gr %>% segment()

is_nucleosome <- rowSums(mcols(gr)$state[, 31:33] == 1) == 3
is_nfr <- rowSums(mcols(gr)$state[, 31:33] == -1) == 3

i <- mcols(gr)$state[, 32] == 1
G_max <- 7
clust <- kmeans(mcols(gr)$z_score, G_max)$cluster
lapply(1:G_max, function(i) plot(colMeans(mcols(gr)$z_score[clust == i, ])))

w <- width(gr)[1]
s <- getSeq(Mmusculus, gr)
dinuc <- expand.grid(c('A', 'C', 'G', 'T'), c('A', 'C', 'G', 'T'))
dinuc <- sprintf('%s%s', dinuc[, 1], dinuc[, 2])


par(mfcol = c(5, 7))
for (h in 1:G_max){
Y <- do.call('cbind', lapply(1:(w - 1), function(i){
	table(factor(as.character(subseq(s[clust == h, ], start = i, width = 2)), dinuc)) / sum(clust == h)
}))
image(matrix(colMeans(as.matrix(mcols(gr)$counts[clust == h, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = h)
image(matrix(colMeans(as.matrix(mcols(gr)$fitted_counts[clust == h, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = h)
plot(colMeans(mcols(gr)$coverage[clust == h,]), type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(colSums(Y[c('GG', 'GC', 'CG', 'CC'), ]), type = 'b')
plot(colSums(Y[c('AA', 'AT', 'TA', 'TT'), ]), type = 'b')
}


is_nucleosome <- rowSums(mcols(gr)$z_score[, 32, drop = FALSE] > 0) == 1
is_nfr <- rowSums(mcols(gr)$z_score[, 32, drop = FALSE] < -0) == 1
s <- rep('unknown', length(gr))
s[is_nucleosome] <- 'nucleosome'
s[is_nfr] <- 'nfr'
table(s[mcols(gr)$group == 1], s[mcols(gr)$group == 2])
plot(colMeans(mcols(gr)$z_score[is_nucleosome, ]))
plot(colMeans(mcols(gr)$z_score[is_nfr, ]))
image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nucleosome, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
image(matrix(colMeans(as.matrix(mcols(gr)$counts[is_nfr, ])), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)



# -----------------------------------------------------------------------------------
# [2019-08-07] Read the chemically mapped nucleosome poisitions for mESC
# -----------------------------------------------------------------------------------
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows)

devtools::load_all('packages/compbio')
nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))	# a simple seqnames/start/end format
gr_nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
gr_nuc <- GRanges(seqnames = gr_nuc[, 1], range = IRanges(gr_nuc[, 2], gr_nuc[, 3]))
gr_nuc <- add.seqinfo(gr_nuc, 'mm10')


#ncp <- rtracklayer::import(sprintf('%s/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz', sra.run.dir('GSM2183909')), format = 'BED', which = reduce(windows))
#ncp <- resize(ncp, width = 1, fix = 'center')
#values(ncp)$name <- as.numeric(values(ncp)$name)
#gr_nuc <- ncp[values(ncp)$name > quantile(values(ncp)$name, 0.90)]

#win <- resize(gr, fix = 'center', width = 100); win <- add.seqinfo(win , 'mm10')
win <- gr; win <- add.seqinfo(win , 'mm10')

# --- read the nucleosome predicted by NucleoATAC
dataset <- 'dataset=Maza_version=20170302a'
base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
fs <- sprintf('%s.nucleoatac_%s.occpeaks.bed.gz', base.name, c('chr1', 'chr2', 'chr3'))
gr1 <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
#########gr1 <- gr1[gr1[, 4] > 0.95, ]
gr1 <- GRanges(seqnames = gr1[, 1], range = IRanges(gr1[, 2], gr1[, 3]))
gr1 <- gr1[gr1 %over% win]
gr1 <- add.seqinfo(gr1, 'mm10')

X <- mcols(gr)$fitted_counts; dim(X) <- c(nrow(X), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals)
breaks <- seq(metadata(gr)$fragment_size_range[1], metadata(gr)$fragment_size_range[2], by = metadata(gr)$fragment_size_interval)
is_nfr <- breaks[-1] > 0 & breaks[-1] <= 100
is_mono_nucleosome <- breaks[-1] >= 180 & breaks[-1] <= 247
is_di_nucleosome <- breaks[-1] >= 315 & breaks[-1] <= 473
is_tri_nucleosome <- breaks[-1] >= 558 & breaks[-1] <= 615
mcols(gr)$nucleosome <- rowMax(rowMeans(X[, , is_mono_nucleosome | is_di_nucleosome | is_tri_nucleosome], dims = 2))
R <- log(rowMeans(rowMeans(X[, , is_mono_nucleosome | is_di_nucleosome | is_tri_nucleosome], dims = 2)) + 0.1) - log(rowMeans(rowMeans(X[, , is_nfr], dims = 2)) + 0.1)
gr2 <- gr[mcols(gr)$nucleosome > quantile(mcols(gr)$nucleosome, 0.6)]
#gr2 <- gr[R > quantile(R, 0.8)]


source('analysis/seatac/helper.r'); evaluate_nucleosome_prediction(win, resize(gr2, fix = 'center', width = 1), gr_nuc[gr_nuc %over% win])
source('analysis/seatac/helper.r'); evaluate_nucleosome_prediction(win, resize(gr1, fix = 'center', width = 1), gr_nuc[gr_nuc %over% win])



ncp <- rtracklayer::import(sprintf('%s/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz', sra.run.dir('GSM2183909')), format = 'BED', which = reduce(win))
ncp <- add.seqinfo(ncp, 'mm10')
ncp <- resize(ncp, width = 1, fix = 'center')
values(ncp)$name <- as.numeric(values(ncp)$name)
cvg <- coverage(ncp, weight = 'name')

NCP <- mean(cvg[unlist(slidingWindows(win, width = bin_size, step = bin_size))])
NCP <- matrix(NCP, nrow = length(win), ncol =  metadata(win)$n_bins_per_window, byrow = TRUE)

X1 <- mean(coverage(gr1)[trim(unlist(slidingWindows(win, width = bin_size, step = bin_size)))])
X1 <- matrix(X1, nrow = length(win), ncol =  metadata(win)$n_bins_per_window, byrow = TRUE)

bw_file <- '/panfs/roc/scratch/gongx030/seatac/MNase_mESC.bw'
cvg <- rtracklayer::import(bw_file, which = trim(reduce(windows)), as = 'RleList')
X3 <- mean(cvg[trim(unlist(slidingWindows(windows, width = bin_size, step = bin_size)))])
X3 <- matrix(X3, nrow = length(windows), ncol =  metadata(windows)$n_bins_per_window, byrow = TRUE)

par(mfcol = c(9, 6))
i <- 2
yy <- c(50, 200, 400, 600, 670)
image(matrix(as.matrix(mcols(gr)$counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy, las = 2)
axis(1, at = c(0, 0.5, 1))
image(matrix(as.matrix(mcols(gr)$fitted_counts[i, ]), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals), col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
axis(2, at = (yy - 50) / (670 - 50), label = yy)
axis(1, at = c(0, 0.5, 1))
plot(mcols(gr)$nucleosome[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i", main = 'nucleosome')
plot(mcols(gr)$coverage[i,], type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(mcols(gr)$fitted_coverage[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i")
plot(mcols(gr)$z_score[i, ], type = 'b', lwd = 2, xaxs="i", yaxs="i"); abline(h = 0)
plot(NCP[i, ], type = 'b', xaxs="i", yaxs="i", main = 'NCP score', xpd = TRUE)
plot(X1[i, ], xaxs="i", yaxs="i", main = 'nucleoatac', xpd = TRUE)
plot(X3[i, ], main = 'MNase', xaxs="i", yaxs="i", xpd = TRUE)


# -----------------------------------------------------------------------------------
# [2019-08-12] Check the NCP score surrounding the TSS
# -----------------------------------------------------------------------------------
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
devtools::load_all('packages/compbio')
gr <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream = 1000, downstream = 1000)
gr <- add.seqinfo(gr, 'mm10')

source('analysis/seatac/helper.r'); ncp <- read_ncp_mESC(which = gr)
cvg <- coverage(ncp, weight = 'name')
X <- as(as(cvg[gr[strand(gr) == '+']], 'RleViews'), 'matrix')


# -----------------------------------------------------------------------------------
# [2019-08-12] Predict the NCP score surrounng TSS
# [2019-08-16] Predict the NCP score surrounng TSS for MEF
# -----------------------------------------------------------------------------------

# --- load the trained model on mESC
gs <- 'Maza_mESC'; ps <- 'Maza_mESC'; latent_dim <- 2; expand <- 1000; window_size <- 320; bin_size <- 10; mr <- 10; mc <- 5; ns <- 1
source('analysis/seatac/helper.r'); model_dir <- model_dir_name(gs, ps, latent_dim, expand, window_size, mr, mc, bin_size, ns)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)

library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
devtools::load_all('packages/compbio')
pmt <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream = 500, downstream = 500)
pmt <- add.seqinfo(pmt, 'mm10')
gs <- 'MEF_NoDox'; ps <- 'MEF_NoDox'; step_size <- 10
source('analysis/seatac/helper.r'); peaks <- read_peaks(ps)
source('analysis/seatac/helper.r'); ga <- read_bam_files(gs, peaks, genome = BSgenome.Mmusculus.UCSC.mm10)
devtools::load_all('analysis/seatac/packages/seatac'); windows <- read_windows(ga, peaks, window_size = window_size, step_size = step_size, bin_size = bin_size)

window_dim <- length(windows_test)
feature_dim <- metadata(windows_test)$n_intervals
input_dim <- metadata(windows_test)$n_bins_per_window

vplot <- mcols(windows_test)$counts %>%
  as.matrix() %>%
  array_reshape(c(window_dim, input_dim, feature_dim, 1))

cvg <- mcols(windows_test)$coverage %>%
  array_reshape(c(window_dim, input_dim, 1))

mcols(windows_test)$label_pred <- model %>% predict(list(vplot, cvg), batch_size = bs, verbose = 1)


