library(futile.logger)
library(GenomicRanges)

devtools::load_all('packages/compbio')
devtools::load_all('analysis/seatac/packages/seatac')

# touch everything in the project dir
# find /panfs/roc/scratch/gongx030/seatac -type f -exec touch {} +

model_dir_name <- function(dataset, latent_dim, window_size, bin_size){
	f <- sprintf('analysis/seatac/models/dataset=%s_latent_dim=%d_window_size=%d_bin_size=%d', dataset, latent_dim, window_size, bin_size)
	flog.info(sprintf('model dir: %s', f))
	f
}

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


#' prepare_windows
#'
prepare_windows <- function(gs, window_size = 320, bin_size = 10, genome){

	if (genome == 'hg19'){
		require(BSgenome.Hsapiens.UCSC.hg19)
		genome2 <- BSgenome.Hsapiens.UCSC.hg19
	}else if (genome == 'mm10'){
		require(BSgenome.Mmusculus.UCSC.mm10)
		genome2 <- BSgenome.Mmusculus.UCSC.mm10
  }

	if (gs == 'MEF_NoDox'){
		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_summits.bed'
		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'
	}else if (gs == 'D1_Dox_Etv2_on_MEF'){
		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'
	}else if (gs == 'D1_Dox_Etv2_on_D0+D1_MEF'){
		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam'
		)
	}else
		stop(sprintf('unknown gs: %s', gs))

	peaks <- macs2.read_summits(bed_file)
	flog.info('removing peaks from chrM')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	peaks <- resize(peaks, fix = 'center', width = window_size)
	peaks <- add.seqinfo(peaks, genome)
	ga <- read_bam(bam_files, peaks, genome = genome2, expand = 2000)
	peaks <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size)

	peaks
} # prepare_windows


