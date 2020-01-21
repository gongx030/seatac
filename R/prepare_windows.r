#' prepare_windows
#'
prepare_windows <- function(
	x,
	window_size = 320, 
	bin_size = 5, 
	bigwig_files = NULL,
	bam_files = NULL,
	blacklist_file = NULL,
	fragment_size_range = c(50, 370),
	fragment_size_interval = 5, 
	genome
){

	is_chrM <- seqnames(x) == 'chrM'
	if (any(is_chrM)){
		flog.info(sprintf('removing %d x on chrM', sum(is_chrM)))
		x <- x[!is_chrM]
	}

	flog.info(sprintf('resizing each peak to %d bp', window_size))
	x <- resize(x, fix = 'center', width = window_size)

	seqlevels(x) <- seqlevels(genome)
	seqlengths(x) <- seqlengths(genome)
	seqinfo(x) <- seqinfo(genome)

	if (!is.null(blacklist_file)){
		if (!file.exists(blacklist_file)){
			flog.error(sprintf('%s does not exist', blacklist_file))
		}else{
			blacklist <- read.table(gzfile(blacklist_file), sep = '\t')
			blacklist <- GRanges(seqnames = blacklist[, 1], range = IRanges(blacklist[, 2], blacklist[, 3]))
			i <- x %over% blacklist
			flog.info(sprintf('removing %d overlaping with the blacklist', sum(i)))
			x <- x[!i]
		}
	}

	if (!is.null(bigwig_files)){
		for (i in 1:length(bigwig_files)){
			if (!file.exists(bigwig_files[i])){
				flog.error(sprintf('%s does not exist', bigwig_files[i]))
			}else{
				flog.info(sprintf('reading %s', bigwig_files[i]))
				gr <- rtracklayer::import(bigwig_files[i], format = 'BigWig', which = reduce(x))
				seqlevels(gr) <- seqlevels(genome)
				seqlengths(gr) <- seqlengths(genome)
				seqinfo(gr) <- seqinfo(genome)
				cvg <- coverage(gr, weight = as.numeric(mcols(gr)$score))
				mcols(x)[[names(bigwig_files)[i]]] <- as(as(cvg[x], 'RleViews'), 'matrix')
			}
		}
	}

	ga <- read_bam(bam_files, x, genome = genome, expand = window_size * 2)

	x <- read_fragment_size_matrix(
		ga, 
		x, 
		window_size = window_size, 
		bin_size = bin_size, 
		fragment_size_range = fragment_size_range, 
		fragment_size_interval = fragment_size_interval
	)
	x

} # prepare_windows
