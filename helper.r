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
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20170519a/MNase_treat_pileup.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, 'mm10')
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')

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

#		nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.nucleoatac.nucleoatac_signal.smooth.bedgraph.gz'
#		flog.info(sprintf('reading %s', nucleoatac_file))
#		nuc <- rtracklayer::import(nucleoatac_file, format = 'BED', which = reduce(peaks))
#		nuc <- resize(nuc, width = 1, fix = 'center')
#		nuc <- add.seqinfo(nuc, 'mm10')
#		cvg <- coverage(nuc, weight = as.numeric(mcols(nuc)$name))
#		mcols(peaks)$nucleoatac_signal <- as(as(cvg[peaks], 'RleViews'), 'matrix')

#		peaks_extend <- resize(peaks, fix = 'center', width = window_size + 73 * 2)

#		ncp_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Voong_version=20190805a/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz'
#		flog.info(sprintf('reading %s', ncp_file))
#		ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(peaks_extend))
#		ncp <- resize(ncp, width = 1, fix = 'center')
#		ncp <- add.seqinfo(ncp, 'mm10')

#		cvg <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
#		X <- as(as(cvg[peaks_extend], 'RleViews'), 'matrix')
#		w <- exp(-((-73:73) / 20)^2 / 2)
#		js <- (73 + 1):(73 + window_size)
#		mcols(peaks)$weighted_nucleosome_score <- do.call('cbind', lapply(js, function(j) rowSums(X[, (j - 73):(j + 73)] %*% diag(w))))
#		mcols(peaks)$min_weighted_nucleosome_score <- rowMins(mcols(peaks)$weighted_nucleosome_score)
#		mcols(peaks)$max_weighted_nucleosome_score <- rowMaxs(mcols(peaks)$weighted_nucleosome_score)
#		mcols(peaks)$mean_weighted_nucleosome_score <- rowMeans(mcols(peaks)$weighted_nucleosome_score)
#		mcols(peaks)$weighted_nucleosome_score <- (mcols(peaks)$weighted_nucleosome_score - mcols(peaks)$min_weighted_nucleosome_score) / (mcols(peaks)$max_weighted_nucleosome_score - mcols(peaks)$min_weighted_nucleosome_score)
#		mcols(peaks)$weighted_nucleosome_score[is.na(mcols(peaks)$weighted_nucleosome_score)] <- 0
#
#		mcols(peaks)$nucleosome_score  <- as(as(cvg[peaks], 'RleViews'), 'matrix')
#		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
#		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
#		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
#		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
#		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0

	}else
		stop(sprintf('unknown gs: %s', gs))

	flog.info('removing peaks from chrM')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	peaks <- resize(peaks, fix = 'center', width = window_size)
	peaks <- add.seqinfo(peaks, genome)

	if (genome == 'mm10'){
		blacklist <- read.table(gzfile('/panfs/roc/scratch/gongx030/datasets/datasets=blacklists_version=20190827a/mm10.blacklist.bed.gz'), sep = '\t')
	}

	blacklist <- GRanges(seqnames = blacklist[, 1], range = IRanges(blacklist[, 2], blacklist[, 3]))
	i <- peaks %over% blacklist
	peaks <- peaks[!i]

	flog.info(sprintf('removing %d peaks overlaping with the blacklist', sum(i)))
	mcols(peaks)$sequence <- getSeq(genome2, peaks)

	ga <- read_bam(bam_files, peaks, genome = genome2, expand = 2000)
	peaks <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size, fragment_size_range = fragment_size_range, fragment_size_interval = fragment_size_interval)

	peaks

} # prepare_windows


