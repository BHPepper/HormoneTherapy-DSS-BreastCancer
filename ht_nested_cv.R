library(caret)
library(MLmetrics)

nested_cv = function(cv_k1 = 10, cv_k2 = 10, seed = 1, 
										 best_perf_method = 'mean', 
										 perf_type = 'low', grid = NULL,  
										 inner_perf_f = NULL, outer_perf_f, score_f,
										 model,
										 dat, response) {
	#' Performs a nested cross-validation to determine performance
	#' with a grid-search of parameters. Generalized to work
	#' for either regression or classification on any performance metric.
	#' 
	#' @param cv_k1 - the number of folds in the model performance CV
	#' @param cv_k2 - the number of folds in the grid-search CV
	#' @param seed - a random seed
	#' @param best_perf_method - 'mean' or 'freq'
	#' @param perf_type - 'low' or 'high', is low or high perf best?
	#' @param grid - matrix or data.frame of parameters, one per column
	#' @param model - a function that returns a fit model
	#' @param inner_perf_f - a function that returns a performance number
	#' @param outer_perf_f - a function that returns a performance number
	#' @param score_f - function that returns number for scoring observations
	#' @param dat - list of data.frames
	#' @param response - string indicating column with response
	#' 
	set.seed(seed)
	
	i_folds = createFolds(dat[[1]][,response], k = cv_k1)
	
	outer_perf = matrix()
	best_params = matrix()
	scores = list()
	y = list()
	dat_i_train = list()
	dat_i_test = list()
	
	for(i in 1:cv_k1) {
		i_train = unlist(i_folds[-i])
		i_test = unlist(i_folds[i])
		
		j_folds = createFolds(dat[[1]][,response][i_train], k = cv_k2)
		
		for(k in 1:length(dat)){
		  dat_i_train[[k]] = dat[[k]][i_train,]
		  dat_i_test[[k]] = dat[[k]][i_test,]
		}
		
		if(!is.null(grid)) {
			inner_perf = list()
			
			for(j in 1:cv_k2) {
				j_train = unlist(j_folds[-j])
				j_test = unlist(j_folds[j])
				
				perf_rows = matrix()
				dat_j_train = list()
				dat_j_test = list()
				
				for(k in 1:length(dat)){
				  dat_j_train[[k]] = dat_i_train[[k]][j_train,]
				  dat_j_test[[k]] = dat_i_train[[k]][j_test,]
				}
				
				# Fit Model to Grid-Search Parameters and Evaluate Performance
				for(row in 1:nrow(grid)) {
					fit = model(dat_j_train, response, grid[row,])
					perf_rows[row] = inner_perf_f(fit, dat_j_test, response)
				}
				
				inner_perf[[j]] = perf_rows
				
			}
			inner_perf = as.matrix(do.call(rbind, inner_perf))
			
			if(best_perf_method == 'mean') {
				if(perf_type == 'low') {
					best = which.min(colMeans(inner_perf, na.rm=T))
				} else {
					best = which.max(colMeans(inner_perf, na.rm=T))
				}
			} else {
				if(perf_type == 'low') {
					best = mfv(apply(inner_perf, 1, which.min))[1]
				} else {
					best = mfv(apply(inner_perf, 1, which.max))[1]
				}
			}
			# Fit on Selected Parameters
			fit = model(dat_i_train, response, grid[best,])
			outer_perf[i] = outer_perf_f(fit, dat_i_test, response)
			scores[[i]] = score_f(fit, dat_i_test, response)
			y[[i]] = dat_i_test[[1]][, response]
			best_params[i] = best
		} else {
		  # Fit on Selected Parameters
		  fit = model(dat_i_train, response)
		  outer_perf[i] = outer_perf_f(fit, dat_i_test, response)
		  scores[[i]] = score_f(fit, dat_i_test, response)
		  y[[i]] = dat_i_test[[1]][, response]
			best_params[i] = NA
		}
	}
	
	return(list(Performance = data.frame(Perf = outer_perf, BestParams = best_params),
							Scores = unlist(scores),
							Y = unlist(y)))
}

