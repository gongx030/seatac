setMethod(
	'save_results',
	signature(
		x = 'VplotsFitted',
		dir = 'character'
	),
	function(
		x,
		dir,
		...
	){

		if (!file.exists(dir))
			dir.create(dir, recursive = TRUE)

		results_file <- sprintf('%s/results.rds', dir)
		flog.info(sprintf('writing %s', results_file))
		saveRDS(x, results_file)

		save_model(x@model, dir)
		
	}
)

