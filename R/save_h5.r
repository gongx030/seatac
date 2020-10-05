setMethod(
	'save_h5',
	signature(
		x = 'VplotsKmers',
		filename = 'character'
	),
	function(x, filename, replace = FALSE){

		if (file.exists(filename)){
			if (!replace){
				stop(sprintf('%s exists.', filename))
			}else
				unlink(filename)
		}

		res <- h5createFile(filename)

		h5createGroup(filename, 'rowData')
		cn <- colnames(rowData(x))
		for (h in cn){
			group <- sprintf('rowData/%s', h)
			y <- rowData(x)[[h]]
			h5write(y, filename, group)
		}

		# save coordinates
		h5createGroup(filename, 'rowRanges')

		h5write(as.character(seqnames(x)), filename, 'rowRanges/seqnames')
		h5write(start(x), filename, 'rowRanges/start')
		h5write(end(x), filename, 'rowRanges/end')

		# save assays
		h5createGroup(filename, 'assays')
		cn <- names(assays(x))
		for (h in cn){
			group <- sprintf('assays/%s', h)
			y <- assays(x)[[h]]
			if (class(y) == 'dgCMatrix'){
				y <- as.matrix(summary(y))
				class(y) <- 'integer'
			}
			h5write(y, filename, group)
		}

		h5createGroup(filename, 'slots')
		cn <- slotNames(x)
		cn <- cn[!cn %in% c('rowRanges', 'colData', 'assays', 'NAMES', 'elementMetadata', 'metadata')]
		for (h in cn){
			group <- sprintf('slots/%s', h)
			h5write(slot(x, h), filename, group)
		}
	}
)

