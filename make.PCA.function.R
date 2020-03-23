#' Perform Principle Component Analysis
#'
#' @param gen.data A numeric genotype matrix (rows = taxa, columns = markers).
#' @return A large prcomp ($sdev, $rotation, $center, $scale, $x)

library(tidyverse)
make.pca <- function(gen.data){
	if (class(gen.data[,1]) != "integer" & class(gen.data[,1]) != "numeric"){
		gen.data <- column_to_rownames(gen.data, var = colnames(gen.data[1]))
	}
	PCA <- prcomp(gen.data)
	return(PCA)
}
