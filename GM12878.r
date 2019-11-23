
# benchmarking cvae on GM12878

library(tensorflow)
tfe_enable_eager_execution(device_policy = 'silent')
library(keras)
library(tfprobability)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)


# -----------------------------------------------------------------------------
# [2019-11-14] Conver the bigwig for GM12878 MNasee data to Wig format
# -----------------------------------------------------------------------------
mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw'
bedgraph_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bedgraph'
command <- sprintf('bigWigToWig %s %s', mnase_file, bedgraph_file)
system(command)

wig_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.wig'
command <- sprintf('perl $HOME/src/bioinfo-tools/bedgraph_to_wig.pl --bedgraph %s --wig %s --step 10', bedgraph_file, wig_file)


# -----------------------------------------------------------------------------
# [2019-11-14] Calling nucleosome peaks using DANPOS on GM12878 MNase-seq data
# https://sites.google.com/site/danposdoc/tutorial/dpos
# -----------------------------------------------------------------------------
wig_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.wig'
output_dir <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Siag_danpos'
command <- sprintf('danpos dpos %s -o %s', wig_file, output_dir)
system(command)


# -----------------------------------------------------------------------------
# [2019-11-14] Training SEATAC model for GM12878
# -----------------------------------------------------------------------------
gs_train <- 'GM12878'; flanking <- 1000; bin_size <- 5; fs <- c(50, 370, 5); genome <- 'hg19'; genome2 <- BSgenome.Hsapiens.UCSC.hg19

# --- prepare windows
source('analysis/seatac/prepare_windows.r'); windows <- prepare_windows(gs_train, flanking, bin_size, genome = genome, fragment_size_range = fs[1:2], fragment_size_interval = fs[3])
source('analysis/seatac/windows.r'); save_windows(windows, gs_train, flanking, bin_size, fs[1:2], fs[3])

# --- load training windows
source('analysis/seatac/windows.r'); windows <- load_windows(gs_train, flanking, bin_size, fs[1:2], fs[3])
#latent_dim <- 10; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 5; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 15; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 20; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 25; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 640; epochs <- 100; bs <- 64; spe <- 20
latent_dim <- 10; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 320; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 5 ; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 160; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 10 ; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 160; epochs <- 100; bs <- 64; spe <- 20
#latent_dim <- 15 ; sequence_dim <- 32; n_components <- 15; mr <- 20; window_size <- 160; epochs <- 100; bs <- 64; spe <- 20

devtools::load_all('analysis/seatac/packages/seatac'); model <- seatac(windows, latent_dim = latent_dim, window_size = window_size, min_reads_per_window = mr, epochs = epochs, batch_size = bs, steps_per_epoch = spe, sequence_dim = sequence_dim, n_components = n_components, genome = genome2)
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs_train, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, n_components, mr, flanking, epochs, bs, spe)
devtools::load_all('analysis/seatac/packages/seatac'); save_model(model, model_dir)

# -----------------------------------------------------------------------------
# [2019-11-14] Load the per-trained model
# -----------------------------------------------------------------------------
source('analysis/seatac/model_dir_name.r'); model_dir <- model_dir_name(gs_train, latent_dim, window_size, bin_size, fs[1:2], fs[3], sequence_dim, n_components, mr, flanking, epochs, bs, spe)
devtools::load_all('analysis/seatac/packages/seatac'); model <- load_model(model_dir)


core_width <- 50
core_bin <- round(metadata(model$train_data)$n_bins_per_window / 2 - (core_width / 2) / bin_size):round(metadata(model$train_data)$n_bins_per_window / 2 + (core_width / 2) / bin_size)
nfr <- rowMeans(model$train_data$nfr[, core_bin])
mono <- rowMeans(model$train_data$mono_nucleosome[, core_bin])
r <- log(mono + 1) - log(nfr + 1)

gr <- resize(granges(model$train_data), fix = 'center', width = core_width)
mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw'
mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(gr))
mnase <- add.seqinfo(mnase, genome)
cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
mnase_seq <- mean(cvg[gr])


# --- MNase-seq profile of nucleosome positions by DANPOS
devtools::load_all('packages/compbio');
nuc_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Siag_danpos/pooled/panfs_roc_scratch_gongx030_datasets_dataset=Kundajie_version=20190802a_GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.smooth.positions.xls'
nuc <- read.table(nuc_file, header = TRUE, sep = '\t')
nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
nuc <- resize(nuc, fix = 'center', width = window_size)
nuc <- add.seqinfo(nuc, genome)
nuc <- trim(nuc)
nuc <- nuc[width(nuc) == window_size]
nuc <- nuc[nuc %over% model$input_data]

mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw'
mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = trim(reduce(nuc)))
mnase <- add.seqinfo(mnase, genome)
cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
M <- as(as(cvg[trim(nuc)], 'RleViews'), 'matrix')


# --- prepare the windows for SeATAC predicting
bam_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_50k_cells.bam'
devtools::load_all('analysis/seatac/packages/seatac'); windows <- make_windows(nuc, bam_file, window_size = window_size, bin_size = bin_size, fragment_size_range = fs[1:2], fragment_size_interval = fs[3], genome = genome2)
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(windows, batch_size = 32)

i <- gr$num_reads > 0
plot(colMeans(log(gr[i]$mono_nucleosome + 1e-5) - log(gr[i]$nfr + 1e-5)), main = model$latent_dim)

plot(colMeans(log(gr$mono_nucleosome + 1e-5)) - colMeans(log(model$train_data$mono_nucleosome + 1e-5)))
plot(colMeans(gr$nfr) - colMeans(model$train_data$nfr))

plot(colMeans(log(gr$mono_nucleosome + 1e-5) - colMeans(log(model$train_data$mono_nucleosome + 1e-5))))
plot(colMeans(log(gr$nfr + 1e-5) - colMeans(log(model$train_data$nfr + 1e-5))))

r <- log(gr$mono_nucleosome + 1e-5) - colMeans(log(model$train_data$mono_nucleosome + 1e-5)) - (log(gr$nfr + 1e-5) - colMeans(log(model$train_data$nfr + 1e-5)))
plot(colMeans(r))


plot(colMeans(gr$nfr))
plot(colMeans(log(gr$mono_nucleosome + 1e-5) - log(gr$nfr + 1e-5)))


# --- Look at tdata in NucleoATAC
nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_20170314b.nucleoatac_signal.smooth.bedgraph.gz'
nucleoatac <- rtracklayer::import(nucleoatac_file, format = 'BED', which = reduce(nuc))
nucleoatac <- resize(nucleoatac, width = 1, fix = 'center')
nucleoatac <- add.seqinfo(nucleoatac, genome)
cvg <- coverage(nucleoatac, weight = as.numeric(mcols(nucleoatac)$name))
M2 <- as(as(cvg[nuc], 'RleViews'), 'matrix')

