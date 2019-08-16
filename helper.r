library(futile.logger)
library(GenomicRanges)

devtools::load_all('packages/compbio')
devtools::load_all('analysis/seatac/packages/seatac')

PROJECT_DIR <- '/panfs/roc/scratch/gongx030/seatac'

# touch everything in the project dir
# find /panfs/roc/scratch/gongx030/seatac -type f -exec touch {} +

list_bam_files <- function(){

	bam <- c(
		'MEF_NoDox' = 'ATAC_MEF_NoDox.bam',
		'MEF_Dox_D1' = 'ATAC_MEF_Dox_D1.bam',
		'MEF_Dox_D2' = 'ATAC_MEF_Dox_D2.bam',
		'MEF_Dox_D7' = 'ATAC_MEF_Dox_D7.bam',
		'MEF_Dox_D7_Flk1pos' = 'ATAC_MEF_Dox_D7_Flk1pos.bam',
		'EB_NoDox_D25' = 'ATAC_EB_NoDox_D25.bam',
		'EB_Dox_D25_Flk1pos' = 'ATAC_EB_Dox_D25_Flk1pos.bam',
		'EB_Dox_D25' = 'ATAC_EB_Dox_D25.bam',
		'Maza_MEF' = 'Maza_MEF.bam'
	)
	filenames <- sprintf('analysis/seatac/data/%s', bam)
	names(filenames) <- names(bam)

	filenames <- c(filenames,
		'Maza_mESC' = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.bam'
	)

	if (!all(file.exists(filenames)))
		stop('BAM files missing')

	filenames

} # list_bam_files

list_peak_files <- function(){

	peaks <- c(
		'MEF_NoDox' = 'ATAC_MEF_NoDox_summits.bed',
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

	if (!all(file.exists(filenames)))
		stop('Peak files missing')

	filenames

}

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

	}else if (ps == 'MEF_active_TSS2'){

		library(TxDb.Mmusculus.UCSC.mm10.knownGene)
		x <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, upstream = 1000, downstream = 1000)
		x <- granges(x)
		names(x) <- NULL
		x <- unique(x)

	}else
		stop(sprintf('unknown ps: %s', ps))
	x
		
}


read_bam_files <- function(gs, peaks, genome, expand = 2000){

	ga_file <- sprintf('%s/data/summits/%s_expand=%d_ga.rds', PROJECT_DIR, paste(gs, collapse = '+'), expand)
	if (!file.exists(ga_file)){
		x <- readBAM(list_bam_files()[gs], peaks, genome, expand = expand)
		flog.info(sprintf('writing %s', ga_file))
		saveRDS(x, ga_file)
	}	
	flog.info(sprintf('reading %s', ga_file))
	x <- readRDS(ga_file)
	x

} # read_bam_files


read_windows <- function(ga, peaks, window_size = 320, step_size = 32, bin_size = 10, exclude_exons = FALSE, txdb, min_reads_per_window = 15){

	peaks <- resize(peaks, fix = 'center', width = 1)

	devtools::load_all('analysis/seatac/packages/seatac')
	flog.info(sprintf('# total input peaks: %d', length(peaks)))

	flog.info('removing chrM reads')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	peaks <- add.seqinfo(peaks, 'mm10')

	peaks <- resize(peaks, fix = 'center', width = window_size)
	peaks <- unlist(slidingWindows(peaks, step = step_size, width = step_size))
	peaks <- resize(peaks, fix = 'center', width = window_size)

	gr <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size)
	
	flog.info(sprintf('# total windows: %d', length(gr)))

	A <- matrix(mcols(gr)$num_reads >= min_reads_per_window, nrow = length(peaks), ncol = metadata(ga)$num_samples)
	i <- rowSums(A) == metadata(ga)$num_samples
	gr <- gr[i]
	flog.info(sprintf('# total windows more than %d PE reads in all %d samples: %d', min_reads_per_window, metadata(ga)$num_samples, sum(i)))

	if (exclude_exons){
		flog.info('removing windows overlapping with exons')
		ex <- exons(txdb)
		is_exon <- 1:length(gr) %in% as.matrix(findOverlaps(gr, ex))[, 1]
		gr <- gr[!is_exon]
		flog.info(sprintf('# total windows: %d', length(gr)))
	}
	gr

} # read_windows


get_seatac_res <- function(gs, window_size, bin_size, latent_dim, n_components, prior, beta, epochs, batch_effect){

	output_dir <- sprintf('%s/data/summits/%s_window=%d_bin=%d', PROJECT_DIR, paste(gs, collapse = '+'), window_size, bin_size)
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

	output_dir <- sprintf('%s/data/summits/%s_window=%d_bin=%d', PROJECT_DIR, paste(gs, collapse = '+'), window_size, bin_size)
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

model_dir_name <- function(dataset, peakset, latent_dim, expand, window_size, min_reads_per_window, min_reads_coverage, bin_size, negative_sample_ratio){
	 f <- sprintf('analysis/seatac/models/dataset=%s_peakset=%s_latent_dim=%d_expand=%d_window_size=%d_min_reads_per_window=%d_min_reads_coverage=%d_bin_size=%d_negative_sample_ratio=%d', dataset, peakset, latent_dim, expand, window_size, min_reads_per_window, min_reads_coverage, bin_size, negative_sample_ratio)
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


prepare_training_windows <- function(peaks, expand = 2000, step_size = 10, window_size = 320, negative_sample_ratio = 1){

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

	# the NFR should be any region that are outside of center 147 bp of any known nuc center
	nfr <- setdiff(peaks, resize(nuc, fix = 'center', width = 147))
	nfr <- unlist(slidingWindows(nfr, width = window_size, step = step_size))
	nfr <- nfr[width(nfr) == window_size]
	nfr <- nfr[nfr %within% peaks]

	# reading the NCP score on the candidate NFR region
	ncp <- read_ncp_mESC(which = reduce(peaks))
	cvg <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
	nfr <- nfr[resize(nfr, width = 147, fix = 'center') %over% ncp]
	mcols(nfr)$ncp_score <- mean(cvg[resize(nfr, width = 147, fix = 'center')])
	nfr <- nfr[order(mcols(nfr)$ncp_score)[1:n_neg]]

	gr <- c(nuc, nfr)
	mcols(gr)$ncp_score <- mean(cvg[resize(gr, width = 147, fix = 'center')])
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


find_nfr_mESC <- function(){
	ncp <- read_ncp_mESC()
	browser()
}
