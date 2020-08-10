#' 
#'
setMethod(
	'save_hdf5',
	signature(
		x = 'Vplots',
		file = 'character'
	),
  function(x, file, ...){

		if (file.exists(file))
			unlink(file)

		h5createFile(file)


		y <- granges(x) %>%
			as.data.frame()

#		h5createGroup(file, 'ranges')
		h5write(y, file, 'ranges')

		mcols(x) %>%
			as.data.frame()
		h5createFile(file)
		h5write(mcols(x), file, 'mcols/counts')

		browser()



		slotNames(x)



#		 [1] "kmers"                  "k"                      "fragment_size_range"    "fragment_size_interval" "bin_size"               "window_size"            "n_intervals"            "n_bins_per_window"      "breaks"                 "centers"                "positions"              "seqnames"
#		[13] "ranges"                 "strand"                 "seqinfo"                "elementMetadata"        "elementType"            "metadata"


	}
)

setMethod(
	'save_hdf5',
	signature(
		x = 'VplotsKmers',
		file = 'character'
	),
  function(x, file, ...){
		callNextMethod()
		print('there')
	}
)

