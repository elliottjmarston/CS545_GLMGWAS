#' GLM Function
#' 
#' @param geno A matrix of genotype data with dimensions n x m 
#' @param pheno A matrix of phenotype data with dimensions n x 1
#' @param covariates A matrix of covariate data with dimensions n x t

GLM.func <- function(geno = NULL, pheno = NULL, covariates = NULL){
  working.geno <- geno.data[,-1] #pull user input into genotype matrix, possible data transformation?
  working.pheno <- pheno.data[,-1] #pull user input into phenotype matrix
  working.cov <- covariate.data[,-1] #pull user input into covariate matrix
  
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
      coeff.matrix <- as.numeric(unlist(cbind(1, working.cov, SNP))) #creates matrix of X values to calculate b on
      coeff.matrix <- as.matrix(coeff.matrix)
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
} #end function

#' additional notes: may end up using this as the base of our major function 
#' and call our other functions above the beginning of the for-loop?
