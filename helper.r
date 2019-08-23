library(futile.logger)
library(GenomicRanges)

devtools::load_all('packages/compbio')
devtools::load_all('analysis/seatac/packages/seatac')

# touch everything in the project dir
# find /panfs/roc/scratch/gongx030/seatac -type f -exec touch {} +

#' list_peak_files
#'
list_peak_files <- function(){

	peaks <- c(
		'MEF_Dox_D1' = 'ATAC_MEF_Dox_D1_summits.bed',
		'MEF_Dox_D2' = 'ATAC_MEF_Dox_D2_summits.bed',
		'MEF_Dox_D7' = 'ATAC_MEF_Dox_D7_summits.bed',
		'MEF_Dox_D7_Flk1pos' = 'ATAC_MEF_Dox_D7_Flk1pos_summits.bed',
		'EB_NoDox_D25' = 'ATAC_EB_NoDox_D25_summits.bed',
		'EB_Dox_D25_Flk1pos' = 'ATAC_EB_Dox_D25_Flk1pos_summits.bed',
		'EB_Dox_D25' = 'ATAC_EB_Dox_D25_summits.bed',
		'Maza_MEF' = 'Maza_MEF_summits.bed'
	)
	filenames <- sprintf('analysis/seatac/data/%s', peaks)
	names(filenames) <- names(peaks)

	filenames <- c(filenames, c(
		'MEF_Dox_D1_Etv2' = sprintf('%s/etv2_pioneer/macs2/MEF_Dox_D1_Etv2_summits.bed', Sys.getenv('TMPDIR')),	# Etv2 ChiP-seq
		'Maza_mESC' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed',
		'Maza_mESC_chr1' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed',
		'Maza_mESC_chr2' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed',
		'Maza_mESC_chr3' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed',
		'Maza_mESC_chr1-3' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed',
		'Maza_mESC_chr1-19' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed'
	))

	filenames <- c(filenames, c(
		'MEF_NoDox' = sprintf('%s/MEF_NoDox_summits.bed', dataset_dir('dataset=Etv2ATAC_version=20190228a')),
		'GM12878_50k_cells' = sprintf('%s/GM12878_50k_cells_summits.bed', dataset_dir('dataset=Buenrostro_version=20170721a'))
	))

	if (!all(file.exists(filenames)))
		stop('Peak files missing')

	filenames

}

#' read_peaks
#'
read_peaks <- function(ps){

	peak_files <- list_peak_files()

	if (ps %in% names(peak_files)){

		flog.info(sprintf('reading peak file %s', peak_files[ps]))
		x <- read.table(peak_files[ps], header = FALSE, sep = '\t')
		x <- GRanges(seqnames = x[, 1], range = IRanges(x[, 2], x[, 3]))

		if (ps == 'Maza_mESC_chr2'){
			x <- x[seqnames(x) == 'chr2']
		}else if (ps == 'Maza_mESC_chr1'){
			x <- x[seqnames(x) == 'chr1']
		}else if (ps == 'Maza_mESC_chr3'){
			x <- x[seqnames(x) == 'chr3']
		}else if (ps == 'Maza_mESC_chr1-3'){
			x <- x[seqnames(x) %in% c('chr1', 'chr2', 'chr3')]
		}else if (ps == 'Maza_mESC_chr1-19'){
			x <- x[seqnames(x) %in% sprintf('chr%d', 1:19)]
		}

	}else if (ps %in% c('MEF_active_TSS', 'MEF_active_TSS_0_400', 'MEF_active_TSS_0_800')){

		library(org.Mm.eg.db)
		library(TxDb.Mmusculus.UCSC.mm10.knownGene)

		x <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream = 1000, downstream = 1000)
		dataset <- 'dataset=MEF_RNA-seq_version=20190801a'
		se_file <- sprintf('analysis/datasets/%s.rds', dataset)
		flog.info(sprintf('reading MEF RNA-seq data: %s', se_file))
		se <- readRDS(se_file)
		se <- se[order(assays(se)$counts[, 1], decreasing = TRUE)]
		se <- se[!duplicated(rowData(se)$symbol)]
		bm <- select(org.Mm.eg.db, rowData(se)$symbol, c('ENTREZID'), 'ALIAS')
		x <- x[mcols(x)$gene_id %in% bm[, 'ENTREZID']]
		bm <- bm[bm[, 'ENTREZID'] %in% mcols(x)$gene_id, ]
		se <- se[rowData(se)$symbol %in% bm[, 'ALIAS']]

		A <- sparseMatrix(
			i = as.numeric(factor(bm[, 'ALIAS'], rowData(se)$symbol)),
			j = as.numeric(factor(bm[, 'ENTREZID'], mcols(x)$gene_id)),
			dims = c(length(se), length(x))
		)

		mcols(x)$expression <- t(A) %*% assays(se)$counts[, 1, drop = FALSE]

		if (ps == 'MEF_active_TSS_0_400'){
			x <- narrow(x, start = 1001, end = 1000 + 400)
		}

		if (ps == 'MEF_active_TSS_0_800'){
			x <- narrow(x, start = 1001, end = 1000 + 800)
		}

	}else if (ps == 'mouse_TSS'){

		library(TxDb.Mmusculus.UCSC.mm10.knownGene)
		x <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, upstream = window_size / 2, downstream = window_size / 2)
		x <- granges(x)
		names(x) <- NULL
		x <- unique(x)

	}else if (ps == 'human_TSS'){

		library(TxDb.Hsapiens.UCSC.hg19.knownGene)
		x <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = window_size / 2, downstream = window_size / 2)
		x <- granges(x)

	}else
		stop(sprintf('unknown ps: %s', ps))
	x
		
}


read_bam_files <- function(gs, peaks, genome, expand = 2000){

	if (genome == 'hg19'){
		require(BSgenome.Hsapiens.UCSC.hg19)
		genome2 <- BSgenome.Hsapiens.UCSC.hg19
	}else if (genome == 'mm10'){
		require(BSgenome.Mmusculus.UCSC.mm10)
		genome2 <- BSgenome.Mmusculus.UCSC.mm10
	}

	ga_file <- sprintf('%s/data/summits/%s_expand=%d_ga.rds', project_dir('seatac'), paste(gs, collapse = '+'), expand)
	if (!file.exists(ga_file)){
		x <- read_bam(list_bam_files()[gs], peaks, genome2, expand = expand)
		flog.info(sprintf('writing %s', ga_file))
		saveRDS(x, ga_file)
	}	
	flog.info(sprintf('reading %s', ga_file))
	x <- readRDS(ga_file)
	x

} # read_bam_files


get_seatac_res <- function(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect){

	output_dir <- sprintf('%s/data/summits/%s_window=%d_bin=%d', project_dir('seatac'), paste(gs, collapse = '+'), window_size, bin_size)
	gr_file <- sprintf('%s/latent=%d_components=%d_prior=%s_beta=%d_epochs=%d_batch=%s.rds', output_dir, latent_dim, n_components, prior, beta, epochs, batch_effect)
	if (!file.exists(gr_file))
		stop(sprintf('%s does not exist', gr_file))

	flog.info(sprintf('reading %s', gr_file))
	gr <- readRDS(gr_file)
	metadata(gr)$peak_set <- sprintf('%s_window=%d_bin=%d', paste(gs, collapse = '+'), window_size, bin_size)
	metadata(gr)$gs <- gs
	gr

}


save_seatac_res <- function(gr, gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect){

	output_dir <- sprintf('%s/data/summits/%s_window=%d_bin=%d', project_dir('seatac'), paste(gs, collapse = '+'), window_size, bin_size)
	if (!file.exists(output_dir)) 
		dir.create(output_dir)
	gr_file <- sprintf('%s/latent=%d_components=%d_prior=%s_beta=%d_epochs=%d_batch=%s.rds', output_dir, latent_dim, n_components, prior, beta, epochs, batch_effect)
	if (file.exists(gr_file))
		stop(sprintf('%s exist', gr_file))
	else{
		flog.info(sprintf('writing %s', gr_file))
		saveRDS(gr, file = gr_file)
	}
}

model_dir_name <- function(dataset, latent_dim, window_size, min_reads_per_window, min_reads_coverage, bin_size){
	f <- sprintf('analysis/seatac/models/dataset=%s_latent_dim=%d_window_size=%d_min_reads_per_window=%d_min_reads_coverage=%d_bin_size=%d', dataset, latent_dim, window_size, min_reads_per_window, min_reads_coverage, bin_size)
	flog.info(sprintf('model dir: %s', f))
	f
}

evaluate_nucleosome_prediction <- function(y, y_pred, y_true){

	label <- y %over% y_true
	label_pred <- y %over% y_pred
#	label_pred <- sample(label_pred)
#	label_pred <- rep(TRUE, length(y))

	print(table(label, label_pred))

	tp <- sum(label & label_pred)
	fp <- sum(!label & label_pred)
	fn <- sum(label & !label_pred)
	tn <- sum(!label & !label_pred)
	sensitivity <- tp / (tp + fn)
	specificity <- tn / (tn + fp)

	f1 <- 2 * tp / (2 * tp + fp + fn)
	acc <-(tp + tn) / (tp + tn + fp + fn)
	ppv <- tp / (tp + fp)
	flog.info(sprintf('# False Positive: %d', fp))
	flog.info(sprintf('# True Positive: %d', tp))
	flog.info(sprintf('# True Negative: %d', tn))
	flog.info(sprintf('# False Negative: %d', fn))
	flog.info(sprintf('# TRUE: %d', sum(label)))
	flog.info(sprintf('# FALSE: %d', sum(!label)))
	flog.info(sprintf('sensitivity: %.3f', sensitivity))
	flog.info(sprintf('specificity : %.3f', specificity))
	flog.info(sprintf('PPV: %.3f', ppv))
	flog.info(sprintf('F1: %.3f', f1))
	flog.info(sprintf('Accuracy: %.3f', acc))
}


prepare_training_windows <- function(peaks, expand = 2000, step_size = 10, window_size = 320, core_width = 50, negative_sample_ratio = 1){

	peaks <- resize(peaks, width = expand, fix = 'center')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	peaks <- add.seqinfo(peaks, 'mm10')

	nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))  # a simple seqnames/start/end format
	flog.info(sprintf('reading %s', nc_mm10_gz_file))
	nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
	nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
	nuc <- add.seqinfo(nuc, 'mm10')

	# the 500-bp window centered at the known nucleosome center must be within peak regions
	# so that there will be enough ATAC-seq reads
	nuc <- nuc[nuc %within% peaks]
	nuc <- resize(nuc, width = window_size, fix = 'center')
	nuc <- nuc[nuc %within% peaks]	

	n_pos <- length(nuc)
	n_neg <- n_pos * negative_sample_ratio
	flog.info(sprintf('# positive samples: %d', n_pos))
	flog.info(sprintf('# negative samples: %d', n_neg))
	flog.info(sprintf('# total windows: %d', n_pos + n_neg))

	# the NFR should be any region that are outside of "core_width" region of any known nuc center
	# the core_width was initially set as 147 bp, and thus only use the regions outside as the negative samples
	# However, it appears when using such a trained model for predicting the nucleosome at higher reolustion e.g. 10bp
	# the model tends to over-estimate the # of nucleosome.  I think this might due to that the difference between pos/neg samples
	# are too small
	nfr <- setdiff(peaks, resize(nuc, fix = 'center', width = core_width))
	nfr <- unlist(slidingWindows(nfr, width = window_size, step = step_size))
	nfr <- nfr[width(nfr) == window_size]
	nfr <- nfr[nfr %within% peaks]

	# reading the NCP score on the candidate NFR region
	ncp <- read_ncp_mESC(which = reduce(peaks))
	cvg <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
	nfr <- nfr[resize(nfr, width = core_width, fix = 'center') %over% ncp]
	mcols(nfr)$ncp_score <- mean(cvg[resize(nfr, width = core_width, fix = 'center')])
	nfr <- nfr[order(mcols(nfr)$ncp_score)[1:n_neg]]

	gr <- c(nuc, nfr)
	mcols(gr)$ncp_score <- mean(cvg[resize(gr, width = core_width, fix = 'center')])
	mcols(gr)$label <- rep(c(TRUE, FALSE), c(n_pos, n_neg))
	gr	

} # prepare_training_peaks


read_ncp_mESC <- function(which = NULL){
	if (is.null(which))
		stop('which cannot be NULL')
	ncp_file <- sprintf('%s/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz', sra.run.dir('GSM2183909'))
	flog.info(sprintf('reading %s', ncp_file))
	ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(which))
	ncp <- add.seqinfo(ncp, 'mm10')
	ncp <- resize(ncp, width = 1, fix = 'center')
	values(ncp)$name <- as.numeric(values(ncp)$name)
	ncp
}


read_nucleosome_mESC <- function(){
	nc_mm10_gz_file <- sprintf('%s/GSM2183909_unique.map_95pc_mm10.bed.gz', sra.run.dir('GSM2183909'))  # a simple seqnames/start/end format
	flog.info(sprintf('reading %s', nc_mm10_gz_file))
	nuc <- read.table(gzfile(nc_mm10_gz_file), header = FALSE, sep = '\t')
	nuc <- GRanges(seqnames = nuc[, 1], range = IRanges(nuc[, 2], nuc[, 3]))
	nuc <- add.seqinfo(nuc, 'mm10')
	nuc
} # read_nucleosome_mESC


#' prepare_training_windows2
#'
prepare_training_windows2 <- function(gs, window_size = 320, bin_size = 10, genome){

#	ga <- read_bam_files(gs, peaks, genome = genome)
	if (genome == 'hg19'){
		require(BSgenome.Hsapiens.UCSC.hg19)
		genome2 <- BSgenome.Hsapiens.UCSC.hg19
	}else if (genome == 'mm10'){
		require(BSgenome.Mmusculus.UCSC.mm10)
		genome2 <- BSgenome.Mmusculus.UCSC.mm10
  }

	if (gs == 'Maza_mESC'){

		ncp <- read_ncp_mESC(which = reduce(peaks_extend))
		cvg_ncp <- coverage(ncp, weight = 'name')
		X <- as(as(cvg_ncp[peaks_extend], 'RleViews'), 'matrix')
		w <- exp(-((-73:73) / 20)^2 / 2)
		js <- (73 + 1):(73 + window_size)
		mcols(peaks)$nucleosome_score <- do.call('cbind', lapply(js, function(j) rowSums(X[, (j - 73):(j + 73)] %*% diag(w))))
		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0

		nuc <- read_nucleosome_mESC()
		nuc <- nuc[nuc %over% peaks]
		nuc <- resize(nuc, width = 50, fix = 'center')
		cvg <- coverage(nuc)
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')
		mcols(peaks)$nucleosome_score[mcols(peaks)$nucleosome_score > 1] <- 1

		flog.info('reading NucleoATAC smooth signal')
		dataset <- 'dataset=Maza_version=20170302a'
		base.name <- sprintf('%s/mESC_ATAC', dataset_dir(dataset))
		fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
		na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
		na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
		na <- add.seqinfo(na, 'mm10')
		cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg_nucleoatac[peaks], 'RleViews'), 'matrix')

	}else if (gs == 'MEF_NoDox'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		flog.info('removing peaks from chrM')
		peaks <- peaks[seqnames(peaks) != 'chrM']
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, genome)
		peaks_extend <- resize(peaks, fix = 'center', width = window_size + 73 * 2)

		bam_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'
		ga <- read_bam(bam_file, peaks, genome = genome2, expand = 2000)

		bw_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20170519a/MNase_treat_pileup.bw'
		flog.info(sprintf('reading %s', bw_file))
		cvg <- rtracklayer::import(bw_file, which = reduce(peaks_extend), as = 'Rle')
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')
		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0

#		flog.info('reading NucleoATAC smooth signal')
#		dataset <- 'dataset=Etv2ATAC_version=20190228a'
#		base.name <- sprintf('%s/MEF_NoDox', dataset_dir(dataset))
#		fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
#		na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
#		na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
#		na <- add.seqinfo(na, 'mm10')
#		cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))
#		mcols(peaks)$nucleoatac_signal <- as(as(cvg_nucleoatac[peaks], 'RleViews'), 'matrix')

	}else if (gs == 'MEF_D1_Dox_Etv2'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		flog.info('removing peaks from chrM')
		peaks <- peaks[seqnames(peaks) != 'chrM']
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, genome)
		peaks_extend <- resize(peaks, fix = 'center', width = window_size + 73 * 2)

		bam_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'
		ga <- read_bam(bam_file, peaks, genome = genome2, expand = 2000)

		bw_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20170519a/MNase_treat_pileup.bw'
#		bw_file <- '/home/garrydj/gongx030/rlib/analysis/etv2_pioneer/data/Teif/MEF_Coverage_pileup.bw'
		flog.info(sprintf('reading %s', bw_file))
		cvg <- rtracklayer::import(bw_file, which = reduce(peaks_extend), as = 'Rle')
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')
		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0
	
	}else if (gs == 'GM12878_50k_cells'){

		bw_files <- get_bigwig_files()
		bw_file <- bw_files['GM12878_MNase']
		flog.info(sprintf('reading %s', bw_file))
		cvg <- rtracklayer::import(bw_file, which = reduce(peaks_extend))
		cvg <- add.seqinfo(cvg, genome)
		cvg <- coverage(cvg, weight = mcols(cvg)$score)
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')
		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0

		flog.info('reading NucleoATAC smooth signal')
		dataset <- 'dataset=Buenrostro_version=20170721a'
		base.name <- sprintf('%s/GM12878_50k_cells', dataset_dir(dataset))
		fs <- sprintf('%s.nucleoatac_%s.nucleoatac_signal.smooth.bedgraph.gz', base.name, sprintf('chr%s', c(1:19, 'X', 'Y')))
		na <- do.call('rbind', lapply(fs, function(f) read.table(gzfile(f), header = FALSE, sep = '\t')))
		na <- GRanges(seqnames = na[, 1], range = IRanges(na[, 2], na[, 3]), score = na[, 4])
		na <- add.seqinfo(na, genome)
		cvg_nucleoatac <- coverage(na, weight = as.numeric(mcols(na)$score))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg_nucleoatac[peaks], 'RleViews'), 'matrix')

	}else
		stop(sprintf('unknown ps: %s', ps))

	peaks <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size)

	peaks
}

save_windows <- function(windows, gs, window_size, bin_size){
	f <- sprintf('%s/windows/dataset=%s_window_size=%d_bin_size=%d.rds', project_dir('seatac'), gs, window_size, bin_size)
	flog.info(sprintf('writing %s', f))
	saveRDS(windows, f)
}

load_windows <- function(gs, window_size, bin_size){
	f <- sprintf('%s/windows/dataset=%s_window_size=%d_bin_size=%d.rds', project_dir('seatac'), gs, window_size, bin_size)
	if (!file.exists(f))
		stop(sprintf('%s does not exist', f))
	flog.info(sprintf('reading %s', f))
	windows <- readRDS(f)
	windows
} # load_windows
