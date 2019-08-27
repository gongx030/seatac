library(futile.logger)
library(GenomicRanges)

devtools::load_all('packages/compbio')
devtools::load_all('analysis/seatac/packages/seatac')

# touch everything in the project dir
# find /panfs/roc/scratch/gongx030/seatac -type f -exec touch {} +

read_ncp_mESC <- function(which = NULL){
	if (is.null(which))
		stop('which cannot be NULL')
	ncp_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Voong_version=20190805a/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz'
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



#' prepare_windows
#'
prepare_windows <- function(gs, window_size = 320, bin_size = 10, fragment_size_range, fragment_size_interval, genome){

	if (genome == 'hg19'){
		require(BSgenome.Hsapiens.UCSC.hg19)
		genome2 <- BSgenome.Hsapiens.UCSC.hg19
	}else if (genome == 'mm10'){
		require(BSgenome.Mmusculus.UCSC.mm10)
		genome2 <- BSgenome.Mmusculus.UCSC.mm10
  }

	if (gs == 'MEF_NoDox'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_summits.bed'
		peaks <- macs2.read_summits(bed_file)

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'

	}else if (gs == 'D1_Dox_Etv2_on_MEF'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'

	}else if (gs == 'D1_Dox_Etv2_on_D0+D1_MEF'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam'
		)
	}else if (gs == 'Etv2_MEF_reprogramming'){
		bed_files <- c(
			D1 = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed',
			D2 = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D2_Etv2_summits.bed',
			D7 = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D7_Etv2_summits.bed'
		)
		grs <- lapply(bed_files, function(bed_file){
			gr <- macs2.read_summits(bed_file)
			gr <- gr[mcols(gr)$score > 3] # q-value < 0.001
		})
		grs <- lapply(grs, function(gr) resize(gr, width = 20, fix = 'center'))
		names(grs) <- names(bed_files)
		ol <- findOverlapsOfPeaks(grs[['D1']], grs[['D2']], grs[['D7']])
		peaks <- Reduce('c', ol$peaklist)
		G <- do.call('cbind', lapply(grs, function(gr) peaks %over% gr))
		mcols(peaks)$peak_group <- G

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D2.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos.bam'
		)

	}else if (gs == 'Maza_mESC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed'
		peaks <- macs2.read_summits(bed_file)

		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')
		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.bam'

#		peaks_extend <- resize(peaks, fix = 'center', width = window_size + 73 * 2)
#		ncp <- read_ncp_mESC(which = reduce(peaks_extend))
#		cvg_ncp <- coverage(ncp, weight = 'name')
#		X <- as(as(cvg_ncp[peaks_extend], 'RleViews'), 'matrix')
#		w <- exp(-((-73:73) / 20)^2 / 2)
#		js <- (73 + 1):(73 + window_size)
#		mcols(peaks)$nucleosome_score <- do.call('cbind', lapply(js, function(j) rowSums(X[, (j - 73):(j + 73)] %*% diag(w))))
#		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
#		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
#		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
#		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
#		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0

		ncp_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Voong_version=20190805a/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz'
		flog.info(sprintf('reading %s', ncp_file))
		ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(resize(peaks, width = window_size + 10, fix = 'center')))
		ncp <- resize(ncp, width = 1, fix = 'center')
		ncp <- add.seqinfo(ncp, 'mm10')
		cvg <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
		mcols(peaks)$nucleosome_score  <- as(as(cvg[peaks], 'RleViews'), 'matrix')

	}else
		stop(sprintf('unknown gs: %s', gs))

	flog.info('removing peaks from chrM')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	peaks <- resize(peaks, fix = 'center', width = window_size)
	peaks <- add.seqinfo(peaks, genome)
	ga <- read_bam(bam_files, peaks, genome = genome2, expand = 2000)
	peaks <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size, fragment_size_range = fragment_size_range, fragment_size_interval = fragment_size_interval)

	peaks

} # prepare_windows


