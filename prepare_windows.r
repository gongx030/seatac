devtools::load_all('packages/compbio')
devtools::load_all('analysis/seatac/packages/seatac')

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
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20170519a/MNase_treat_pileup.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, 'mm10')
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_D1_Dox_Etv2_20170314a.nucleoatac_signal.smooth.bedgraph.gz'
		flog.info(sprintf('reading %s', nucleoatac_file))
		nuc <- rtracklayer::import(nucleoatac_file, format = 'BED', which = reduce(peaks))
		nuc <- resize(nuc, width = 1, fix = 'center')
		nuc <- add.seqinfo(nuc, 'mm10')
		cvg <- coverage(nuc, weight = as.numeric(mcols(nuc)$name))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'

	}else if (gs == 'D1_Dox_Etv2_all_ATAC'){

	  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

		# D1 Etv2 ChIP-seq peaks
		peak_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(peak_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, genome = 'mm10')

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D2.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos.bam'
		)

	}else if (gs == 'D1_Etv2_on_MEF_D1_D7Flk1pos_ATAC'){

	  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

		# D1 Etv2 ChIP-seq peaks
		peak_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(peak_file)
		peaks <- add.seqinfo(peaks, genome = 'mm10')
		peaks <- peaks[!peaks %over% exons(TxDb.Mmusculus.UCSC.mm10.knownGene)]
#		peaks <- peaks[!peaks %over% promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)]
		peaks <- resize(peaks, fix = 'center', width = 1)
		flog.info(sprintf('get %d intervals after removing exon region', length(peaks)))

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos.bam'
		)


	}else if (gs == 'D2_Dox_Etv2_on_D0+D2_MEF'){

	  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

		# D2 Etv2
		peak_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D2_Etv2_summits.bed'
		peaks <- macs2.read_summits(peak_file)
		peaks <- add.seqinfo(peaks, genome = 'mm10')

		ex <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene)
		peaks <- peaks[!peaks %over% ex]
		flog.info(sprintf('get %d intervals after removing exon region', length(peaks)))

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D2.bam'
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

	}else if (gs == 'MEF_NoDox+MEF_Dox_D7_Flk1pos+EB_Dox_D25_Flk1pos'){

		library(TxDb.Mmusculus.UCSC.mm10.knownGene)

		bed_files <- c(
			MEF_NoDox = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_summits.bed',
			MEF_Dox_D7_Flk1pos = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos_summits.bed',
			EB_Dox_D25_Flk1pos = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/EB_Dox_D25_Flk1pos_summits.bed'
		)

		grs <- lapply(bed_files, function(bed_file) gr <- macs2.read_summits(bed_file))
		grs <- lapply(grs, function(gr) resize(gr, width = 20, fix = 'center'))
		names(grs) <- names(bed_files)
		ol <- findOverlapsOfPeaks(grs[['MEF_NoDox']], grs[['MEF_Dox_D7_Flk1pos']], grs[['EB_Dox_D25_Flk1pos']])
		peaks <- Reduce('c', ol$peaklist)
		G <- do.call('cbind', lapply(grs, function(gr) peaks %over% gr))
		mcols(peaks)$peak_group <- G

		peaks <- peaks[!peaks %over% exons(TxDb.Mmusculus.UCSC.mm10.knownGene)]
#		peaks <- peaks[!peaks %over% promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)]

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/EB_Dox_D25_Flk1pos.bam'
		)

	}else if (gs == 'D1_Etv2_on_MEF_D7Flk1pos_ATAC'){

	  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

		# D1 Etv2 ChIP-seq peaks
		peak_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(peak_file)
		peaks <- add.seqinfo(peaks, genome = 'mm10')
		peaks <- peaks[!peaks %over% exons(TxDb.Mmusculus.UCSC.mm10.knownGene)]
#		peaks <- peaks[!peaks %over% promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)]
		peaks <- resize(peaks, fix = 'center', width = 1)
		flog.info(sprintf('get %d intervals after removing exon region', length(peaks)))

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos.bam'
		)


	}else if (gs == 'Maza_mESC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed'
		peaks <- macs2.read_summits(bed_file)

		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')
		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.bam'

		nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.nucleoatac.nucleoatac_signal.smooth.bedgraph.gz'
		flog.info(sprintf('reading %s', nucleoatac_file))
		nuc <- rtracklayer::import(nucleoatac_file, format = 'BED', which = reduce(peaks))
		nuc <- resize(nuc, width = 1, fix = 'center')
		nuc <- add.seqinfo(nuc, 'mm10')
		cvg <- coverage(nuc, weight = as.numeric(mcols(nuc)$name))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=mESC_MNase_version=20190925a/MNase_mESC.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, 'mm10')
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score2 <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		peaks_extend <- resize(peaks, fix = 'center', width = window_size + 73 * 2)

		ncp_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Voong_version=20190805a/GSM2183909_Chemical_NCPscore_mm10.sorted_merged.txt.gz'
		flog.info(sprintf('reading %s', ncp_file))
		ncp <- rtracklayer::import(ncp_file, format = 'BED', which = reduce(peaks_extend))
		ncp <- resize(ncp, width = 1, fix = 'center')
		ncp <- add.seqinfo(ncp, 'mm10')

		cvg <- coverage(ncp, weight = as.numeric(mcols(ncp)$name))
		X <- as(as(cvg[peaks_extend], 'RleViews'), 'matrix')
		w <- exp(-((-73:73) / 20)^2 / 2)
		js <- (73 + 1):(73 + window_size)
		mcols(peaks)$weighted_nucleosome_score <- do.call('cbind', lapply(js, function(j) Matrix::rowSums(X[, (j - 73):(j + 73)] %*% diag(w))))
		mcols(peaks)$min_weighted_nucleosome_score <- rowMins(mcols(peaks)$weighted_nucleosome_score)
		mcols(peaks)$max_weighted_nucleosome_score <- rowMaxs(mcols(peaks)$weighted_nucleosome_score)
		mcols(peaks)$mean_weighted_nucleosome_score <- rowMeans(mcols(peaks)$weighted_nucleosome_score)

		mcols(peaks)$weighted_nucleosome_score <- (mcols(peaks)$weighted_nucleosome_score - mcols(peaks)$min_weighted_nucleosome_score) / (mcols(peaks)$max_weighted_nucleosome_score - mcols(peaks)$min_weighted_nucleosome_score)
		mcols(peaks)$weighted_nucleosome_score[is.na(mcols(peaks)$weighted_nucleosome_score)] <- 0

		mcols(peaks)$nucleosome_score  <- as(as(cvg[peaks], 'RleViews'), 'matrix')
		mcols(peaks)$min_nucleosome_score <- rowMins(mcols(peaks)$nucleosome_score)
		mcols(peaks)$max_nucleosome_score <- rowMaxs(mcols(peaks)$nucleosome_score)
		mcols(peaks)$mean_nucleosome_score <- rowMeans(mcols(peaks)$nucleosome_score)
		mcols(peaks)$nucleosome_score <- (mcols(peaks)$nucleosome_score - mcols(peaks)$min_nucleosome_score) / (mcols(peaks)$max_nucleosome_score - mcols(peaks)$min_nucleosome_score)
		mcols(peaks)$nucleosome_score[is.na(mcols(peaks)$nucleosome_score)] <- 0

	}else if (gs == 'MEF_ESC'){

		bed_files <- c(
			MEF = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_summits.bed',
			ESC = '/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC_summits.bed'
		)
		grs <- lapply(bed_files, function(bed_file){
			flog.info(sprintf('reading %s', bed_file))									
			macs2.read_summits(bed_file)
		})
		grs <- lapply(grs, function(gr) resize(gr, width = 20, fix = 'center'))
		names(grs) <- names(bed_files)
		ol <- findOverlapsOfPeaks(grs[['MEF']], grs[['ESC']])
		peaks <- Reduce('c', ol$peaklist)
		G <- do.call('cbind', lapply(grs, function(gr) peaks %over% gr))
		mcols(peaks)$peak_group <- G

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.bam'
		)

	}else if (gs == 'MEF_ESC_with_MEF_Peaks'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_summits.bed'
		peaks <- macs2.read_summits(bed_file)

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Maza_version=20170302a/mESC_ATAC.bam'
		)

	}else if (gs == 'GM12878'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_50k_cells_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, genome)

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, genome)
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_20170314b.nucleoatac_signal.smooth.bedgraph.gz'
		flog.info(sprintf('reading %s', nucleoatac_file))
		nuc <- rtracklayer::import(nucleoatac_file, format = 'BED', which = reduce(peaks))
		nuc <- resize(nuc, width = 1, fix = 'center')
		nuc <- add.seqinfo(nuc, genome)
		cvg <- coverage(nuc, weight = as.numeric(mcols(nuc)$name))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_50k_cells.bam'

	}else if (gs == 'D1_Dox_Etv2_MEF_D1_ATAC'){

		# D1 Etv2 summits on MEF/D1 ATAC-seq

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		bam_files <- c(
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam',
			'/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam'
		)

	}else if (gs == 'D1_Dox_Etv2_D1_ATAC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D1_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1.bam'

	}else if (gs == 'D2_Dox_Etv2_D2_ATAC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D2_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D2.bam'

	}else if (gs == 'D7_Dox_Etv2_D7_ATAC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D7_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7.bam'

	}else if (gs == 'D7_Dox_Etv2_D7Flk1pos_ATAC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ChIPseq_version=20190307a/MEF_Dox_D7_Etv2_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos.bam'

	}else if (gs == 'EB25_Etv2_Flk1pos_ATAC'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=ChoiEtv2ChIP-seq_version=20190909a/etv2_chipseq_polyab_mm10.sorted.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, 'mm10')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/EB_Dox_D25_Flk1pos.bam'

	}else if (gs == 'GM12878_STAT3'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=GM12878_version=20191007a/ENCFF014RBU.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, genome)

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, genome)
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=GM12878_version=20191007a/GM12878_20170314c.nucleoatac_signal.smooth.bw'
		flog.info(sprintf('reading %s', nucleoatac_file))
		nuc <- rtracklayer::import(nucleoatac_file, format = 'BigWig', which = reduce(peaks))
		nuc <- add.seqinfo(nuc, genome)
		cvg <- coverage(nuc, weight = as.numeric(mcols(nuc)$score))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_50k_cells.bam'

	}else if (gs == 'GM12878_ETV6'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=GM12878_version=20191007a/ENCFF692GBM.bed.gz'
		flog.info(sprintf('reading %s', bed_file))
		peaks <- read.table(gzfile(bed_file), header = FALSE, sep = '\t')
		peaks <- GRanges(seqnames = peaks[, 1], range = IRanges(peaks[, 2], peaks[, 3]))
		peaks <- resize(peaks, fix = 'center', width = window_size)
		peaks <- add.seqinfo(peaks, genome)

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Kundajie_version=20190802a/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, genome)
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		nucleoatac_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=GM12878_version=20191007a/GM12878_20170314c.nucleoatac_signal.smooth.bw'
		flog.info(sprintf('reading %s', nucleoatac_file))
		nuc <- rtracklayer::import(nucleoatac_file, format = 'BigWig', which = reduce(peaks))
		nuc <- add.seqinfo(nuc, genome)
		cvg <- coverage(nuc, weight = as.numeric(mcols(nuc)$score))
		mcols(peaks)$nucleoatac_signal <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Buenrostro_version=20170721a/GM12878_50k_cells.bam'

	}else if (gs == 'PHF7_MEF'){

		bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Phf7ChIPseq_version=20190925a/PHF7_summits.bed'
		peaks <- macs2.read_summits(bed_file)
		peaks <- resize(peaks, fix = 'center', width = window_size)

		mnase_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20170519a/MNase_treat_pileup.bw'
		flog.info(sprintf('reading %s', mnase_file))
		mnase <- rtracklayer::import(mnase_file, format = 'BigWig', which = reduce(peaks))
		mnase <- add.seqinfo(mnase, 'mm10')
		cvg <- coverage(mnase, weight = as.numeric(mcols(mnase)$score))
		mcols(peaks)$nucleosome_score <- as(as(cvg[peaks], 'RleViews'), 'matrix')

		bam_files <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox.bam'

	}else
		stop(sprintf('unknown gs: %s', gs))

	flog.info('removing peaks from chrM')
	peaks <- peaks[seqnames(peaks) != 'chrM']
	peaks <- resize(peaks, fix = 'center', width = window_size)
	peaks <- add.seqinfo(peaks, genome)

	if (genome == 'mm10'){
		blacklist_file <- '/panfs/roc/scratch/gongx030/datasets/datasets=blacklists_version=20190827a/mm10.blacklist.bed.gz'
	}else if (genome == 'hg19'){
		blacklist_file <- '/panfs/roc/scratch/gongx030/datasets/datasets=blacklists_version=20190827a/hg19.blacklist.bed.gz'
	}

	blacklist <- read.table(gzfile(blacklist_file), sep = '\t')
	blacklist <- GRanges(seqnames = blacklist[, 1], range = IRanges(blacklist[, 2], blacklist[, 3]))
	i <- peaks %over% blacklist
	peaks <- peaks[!i]

	flog.info(sprintf('removing %d peaks overlaping with the blacklist', sum(i)))
	mcols(peaks)$sequence <- getSeq(genome2, peaks)

	ga <- read_bam(bam_files, peaks, genome = genome2, expand = 2000)
	peaks <- readFragmentSizeMatrix(ga, peaks, window_size = window_size, bin_size = bin_size, fragment_size_range = fragment_size_range, fragment_size_interval = fragment_size_interval)

	peaks

} # prepare_windows


