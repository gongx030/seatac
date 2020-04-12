
setMethod(
	'dim',
	signature(
		x = 'sparse_array'
	),
	function(x){
		x@dims
	}
)
