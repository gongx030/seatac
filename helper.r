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



