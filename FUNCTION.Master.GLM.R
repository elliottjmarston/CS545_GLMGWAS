#' Master GLM Function
#' 
#' @param X A matrix of genotype data with dimensions n x m 
#' @param y A matrix of phenotype data with dimensions n x 1
#' @param C A matrix of covariate data with dimensions n x t
#' @param r Correlation threshold value for filtering principle components colinear to covariates. Automatically set to 0.2
#' @param PC.num the number of principle components to include from the PCA portion of this function's calculation
#' @return A vector of p-values of length m
#' 
#' @examples  
#' \dontrun{
#' Master.GLM.func(X = genotype, y = phenotype, C = covariates)
#' Master.GLM.func(X = genotype, y = phenotype, C = covariates, r = 0.2)
#' }
#' @export
#' 
Master.GLM.func <- function(X = NULL, y = NULL, C = NULL, r = 0.2, PC.num = 10){
  working.geno <- X #pull user input into genotype matrix
  working.pheno <- y[,-1] #pull user input into phenotype matrix
  working.cov <- C #pull user input into covariate matrix
  r.thresh <- r
  PC.n <- PC.num
  
  { #user input formatting block
    if (class(working.geno[,1]) != "integer" & class(working.geno[,1]) != "numeric"){
      working.geno <- column_to_rownames(working.geno, var = colnames(working.geno[1]))
    } #end genotype formatting
    if(class(working.cov[,1]) != "integer" & class((working.cov)[,1]) != "numeric"){
      working.cov <- column_to_rownames(working.cov, var = colnames(working.cov[1]))
    } #end covariate formatting
  }   
  
  { #PCA calculation and filter block
  PCA.values <- prcomp(working.geno) #using prcromp() to perform PCA
  filtered.PCs <- filter.pca(PCA = PCA.values, covs = working.cov, threshold = r.thresh) #creates matrix of PC's filtered for colinearity with covariates    
  working.PCs <- filtered.PCs$x[,1:PC.num] #pulls x number of PCs sorted in order of explaining the most to least data variation
  } #end PCA calcuation and filter block  
  
  sample.n <- nrow(working.geno) #number of samples in data
  marker.n <- ncol(working.geno) #number of markers in data
  
  p.matrix <- matrix(data = NA, nrow = 1, ncol = marker.n) #creation of blank matrix with 1 row and marker.n number of columns
  
  #for loop through genome
  for (i in 1:marker.n) {
    SNP = working.geno[,i] #SNP value set to i column in working.geno
    if(max(SNP)==min(SNP)){ #if monomorphic marker, p = 1
      p.value = 1
    }
    else{
      coeff.matrix <- as.matrix(cbind(1, working.cov, working.PCs, SNP)) #creates matrix of X values to calculate b on
      X.matrix <- t(coeff.matrix)%*%coeff.matrix #creates X^2 multiplication matrix
      inverse.X <- solve(X.matrix) #creates inverse of X^2 matrix
      Y.matrix <- t(coeff.matrix)%*%working.pheno #creates X*Y matrix
      
      b.effect <- inverse.X%*%Y.matrix #multiple two matrices we just created to calculate b
      yb <- coeff.matrix%*%b.effect #intermediary calculation using matrix multiplication
      e <- working.pheno-yb #calculating residual
      
      pheno.n <- length(working.pheno)
      var.e <- sum(e^2)/(pheno.n-1) #calculating variance of residuals
      var.t <- inverse.X*var.e #calculating variance of t-stat
      t.stat <- b.effect/sqrt(diag(var.t)) #calculating t-stat using b and var.t
      
      p.value <- 2*(1-pt(abs(t.stat), pheno.n-2))
    } #monomorphic marker loop
    p.matrix[i] <- p.value[length(p.value)]
  } #SNP loop
  return(p.matrix)
} #end function

