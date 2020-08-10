setMethod(
	'load_results',
	signature(
		dir = 'character'
	),
	function(
		dir,
		...
	){
		if (!file.exists(dir))
			stop(sprintf('%s does not exist dir', dir))

		results_file <- sprintf('%s/results.rds', dir)
		if (!file.exists(results_file))
			stop(sprintf('%s does not exist', results_file))

		flog.info(sprintf('reading %s', results_file))
		x <- readRDS(results_file)

		x@model <- load_model(dir)
		x
	}
)


