#####################################
# FUNCTIONS TO IMPUTE MISSING VALUES 
# IN NUMERICAL MATRIX
#####################################


# get_impute_parameters: estimate parameters of imputation model from training data
# impute_na_gaussian: auxiliary function (should not be directly called by user)
#					  impute missing values based on multivariate normal distribution
# impute_na: impute missing values based either on multivariate normal model 
#			 or on mean vector

get_impute_parameters <- function(x, method = c("gaussian", "mean")) {
	method <- match.arg(method)
	if (method == "mean") {
		mu <- colMeans(x, na.rm = TRUE)
		parameters <- list(mean = mu)
	} else if (method == "gaussian") {
		s <- prelim.norm(x)
		thetahat <- em.norm(s, showits = FALSE, 
			criterion = sqrt(.Machine$double.eps))
		parameters <- getparam.norm(s, thetahat)
	} 
	return(parameters)
}

impute_NA_gaussian <- function(x, parameters) {
	row_miss <- which(rowSums(is.na(x)) > 0)
	if (length(row_miss) == 0) return(x)
	mu <- parameters[["mu"]]
	sigma <- parameters[["sigma"]]
	p <- length(mu)
	shrinkage <- 0.01
	sigma <- (1 - shrinkage) * sigma + diag(shrinkage * m, p)
	for (i in row_miss) {
		obs <- which(!is.na(x[i,])) 
		miss <- (1:p)[-obs]
		if (length(obs) == 0) {
			x[i,] <- mu
		} else if (length(obs) == 1) {
			x[i,miss] <- mu[miss] + sigma[miss,obs] * 
				((x[i,obs] - mu[obs]) / sigma[obs,obs])
		} else {
			x[i,miss] <- mu[miss] + sigma[miss,obs] %*% 
				solve(sigma[obs,obs], x[i,obs] - mu[obs]) 
		}	
	}
	return(x)
}



impute_NA <- function(x, method = c("gaussian", "mean"), parameters) {
	method <- match.arg(method)	
	if (method == "mean") {
		mu <- parameters[["mean"]]
		miss <- which(is.na(x), TRUE)
		x[miss] <- mu[miss[,"col"]]
	} else if (method == "gaussian") {
		x <- impute_NA_gaussian(x, parameters)
	} 
	return(x)
}

