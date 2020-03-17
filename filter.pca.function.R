#' Filter out PCs colinear with other covariates. 
#'
#'@param PCA The output list from the prcomp() function. 
#'@param covs A matrix/data.frame of covariates (rows = individuals, columns = covariates)

library(tidyverse)

filter.pca <- function(PCA, covs){
	
	# ensure the first column is the first covariate
	if (class(covs[,1]) != "integer" & class(covs[,1]) != "numeric"){
		covs <- column_to_rownames(covs, var = colnames(covs[1]))
	}
	
	num.cov <- length(colnames(covs)) #get number of covariates the user input
	
	#function to return boolean index of PCs that are highly correlated with user-input covariates
	filter.pca <- function(PC, covariate){
		if (cor(PC, covariate) >= 0.5){
			return(F)
		}else{
			return(T)
		}
	}
	
	#filter out PCs that are highly correlated with user-input covariates
	for(i in 1:num.cov){
		to.filter <- apply(PCA$x, 2, filter.pca, covariate=covs[i]) 
		PCA$x <- PCA$x[,to.filter]
	}
	
	return(PCA$x)
}