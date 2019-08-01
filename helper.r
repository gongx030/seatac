library(futile.logger)
library(GenomicRanges)

devtools::load_all('packages/compbio')

PROJECT_DIR <- sprintf('%s/seatac', Sys.getenv('TMPDIR'))

flog.info(sprintf('PROJECT_DIR: %s', PROJECT_DIR))
flog.info(sprintf('BAM dir: %s', 'analysis/seatac/data'))

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

	if (!all(file.exists(filenames)))
		stop('BAM files missing')

	names(filenames) <- names(bam)
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

	if (!all(file.exists(filenames)))
		stop('Peak files missing')

	names(filenames) <- names(peaks)
	filenames

}

read_peaks <- function(ps){

	peak_files <- list_peak_files()
	if (ps %in% names(peak_files)){

		flog.info(sprintf('reading peak file %s', peak_files[ps]))
		x <- read.table(peak_files[ps], header = FALSE, sep = '\t')
		x <- GRanges(seqnames = x[, 1], range = IRanges(x[, 2], x[, 3]))

	}else if (ps == 'MEF_active_TSS'){

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
		bm <- bm[bm[, 'ENTREZID'] %in% mcols(peaks)$gene_id, ]
		se <- se[rowData(se)$symbol %in% bm[, 'ALIAS']]

		A <- sparseMatrix(
			i = as.numeric(factor(bm[, 'ALIAS'], rowData(se)$symbol)),
			j = as.numeric(factor(bm[, 'ENTREZID'], mcols(x)$gene_id)),
			dims = c(length(se), length(x))
		)

		mcols(x)$expression <- t(A) %*% assays(se)$counts[, 1, drop = FALSE]

	}else
		stop(sprintf('unknown ps: %s', ps))
	x
		
}


read_bam_files <- function(gs, peaks, genome, expand = 2000){

	ga_file <- sprintf('%s/data/summits/%s_expand=%d_ga.rds', PROJECT_DIR, paste(gs, collapse = '+'), expand)
	if (!file.exists(ga_file)){
		x <- readBAM(list_files()[gs, 'bam'], peaks, genome, expand = expand)
		flog.info(sprintf('writing %s', ga_file))
		saveRDS(x, ga_file)
	}	
	flog.info(sprintf('reading %s', ga_file))
	x <- readRDS(ga_file)
	x

} # read_bam_files


read_windows <- function(gs, ga, peaks, expand, bin_size = 10, exclude_exons = FALSE, txdb, min_reads_per_window = 15){

	peaks <- resize(peaks, fix = 'center', width = 1)

	devtools::load_all('analysis/seatac/packages/seatac')
	flog.info(sprintf('# total input peaks: %d', length(peaks)))

	flog.info('removing chrM reads')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	seqlevels(peaks, pruning.mode = 'coarse') <- seqlevels(ga)

	window_file <- sprintf('%s/data/summits/%s_expand=%d_bin=%d.rds', PROJECT_DIR, paste(ps, collapse = '+'), expand, bin_size)
	if (!file.exists(window_file)){
		gr <- readFragmentSizeMatrix(ga, peaks, window_size = expand, bin_size = bin_size)
		flog.info(sprintf('writing %s', window_file))
		saveRDS(gr, window_file)
	}
	
	flog.info(sprintf('reading %s', window_file))
	gr <- readRDS(window_file)

	flog.info(sprintf('# total windows: %d', length(gr)))

	gr <- gr[mcols(gr)$num_reads >= min_reads_per_window]
	flog.info(sprintf('# total windows after removing windows with less than %d PE reads: %d', min_reads_per_window, length(gr)))

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

model_dir_name <- function(dataset, peakset, expand, latent_dim, n_components, batch_effect, window_size, step_size, min_reads_per_window, min_reads_coverage){
	 f <- sprintf('analysis/seatac/models/dataset=%s_peakset=%s_expand=%d_latent_dim=%d_n_components=%d_batch_effect=%s_window_size=%d_step_size=%d_min_reads_per_window=%d_min_reads_coverage=%d', paste(dataset, collapse = '+'), paste(peakset, collapse = '+'), expand, latent_dim, n_components, batch_effect, window_size, step_size, min_reads_per_window, min_reads_coverage)
	flog.info(sprintf('model dir: %s', f))
	f
}




