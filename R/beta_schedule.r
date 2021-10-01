#' beta_schedule
#' 
#' Schedule for beta during the VAE training
#' 
#' @param method Scheduling method: 'constant', 'monotonic' or 'cyclical_annealing'
#' @param beta0 initial beta value
#' @param epochs Number of total epochs
#' @param ... Additional arguments
#'
#' @export
#'
beta_schedule <- function(method = 'constant', beta0 = 5e-5, epochs, ...){
	params <- list(...)
	if (method == 'constant'){
		beta <- rep(beta0, epochs)
	}else if (method == 'cyclical_annealing'){
		if (is.null(params[['M']])){
			params[['M']] <- 4L
		}
		if (is.null(params[['R']])){
			params[['R']] <- 0.5
		}
		tau <- ((1:epochs - 1) %% floor(epochs / params[['M']])) / floor(epochs / params[['M']])
		beta <- rep(beta0, epochs)
		beta[tau <= params[['R']]] <- tau[tau <= params[['R']]] * beta0 / params[['R']]

	}else if (method == 'monotonic'){

		if (is.null(params[['R']])){
			params[['R']] <- 0.2
		}
		beta <- beta_schedule(method = 'cyclical_annealing', epochs = epochs, M = 1L, R = params[['R']], beta0 = beta0)

	}else{
		stop(sprintf('unknown method: %s', method))
	}
	beta
}

