#'@title COCA
#'@param data Data matrix or data.frame
#'@param rank rank of decomposition
#'@param verbose Should log-likelihood be printed at every step?
#'@param post_svd Should decomposition be post-processed with svd?
#' @description Computes low-rank COCA decomposition of data. 
#' Can be used with missing data. 
#' @details Copula-based decomposition, where all data is treated as
#' continuous. 
#' 
#' For technical details:
#' 
#' Fang Han and Han Liu. Semiparametric principal component analysis. 
#' In NIPS’12: 
#' Proceedings of the 26th Annual Conference on Neural Information 
#' Processing Systems, pages 171–179, 2012.
#' 
#' http://papers.nips.cc/paper/4809-semiparametric-principal-component-analysis.
#' @return 
#' A list with the following fields:
#' 
#'    - \code{A}: Low-rank representation of rows
#'    
#'    - \code{B}: Low-rank representation of columns
#'    
#'    - \code{fittedEsts}: Low-rank representation of full data matrix
#'    
#' Other fields are for internal use. 
#'@examples 
#' test_data = simProblem()$data
#' fit = coca(test_data, rank = 2)
#'@export
coca = function(data, rank, verbose = F, 
                post_svd = T){
  # Returns proportions
  cop_vals = raw2cop(data, "m+1")
  # Transforming to z-value
  cop_vals = qnorm(cop_vals)
  decomp = pca(cop_vals, rank, post_svd = post_svd)
  A = decomp$A
  B = decomp$B
  theta = decomp$fittedEsts
  
  fittedEsts <- matrix(nrow = nrow(data), 
                       ncol = ncol(data))
  colnames(fittedEsts) <- colnames(data)
  for(i in 1:ncol(fittedEsts)){
    this_median_y <- pnorm(theta[,i])
    fittedEsts[,i] = quantile(x = data[,i], 
                              na.rm = T, 
                              type = 1, 
                              probs = this_median_y)
    # Arguemnt type = 1 means quantile returns inverse of EDF
  }
  
  rownames(A) = rownames(data)
  colnames(A) = paste0("PC", 1:rank)
  rownames(B) = colnames(data)
  colnames(B) = paste0("PC", 1:rank)
  dimnames(fittedEsts) = dimnames(data)
  ans =  list(A = A, B = B,
              fittedEsts = fittedEsts, 
              iters = decomp$iters, 
              sigma = sigma)
  return(ans)
}


### PREPROCESSING TOOLS

# Takes in one column of raw data 
# Outputs the copula value for COCA
col2cop = function(vals, denom_method = "m+1"){
  m = sum(!is.na(vals))
  # "m" method is used for XPCA, "m+1" is used for COCA
  # Slightly different ways of calculating Empirical Distribution Function
  if(denom_method == "m"){denom = m}
  else if(denom_method == "m+1"){denom = m+1}
  else{stop("denom_method invalid. Valid options: 'm' or 'm+1'")}
  ans = rank(vals, ties.method = "average", na.last = "keep")/denom
  return(ans)
}

# Takes in full data matrix
# Outputs full L, R matrices 
raw2cop = function(data, denom_method = "m+1"){
  data = as.matrix(data)
  # ans will be the copula values
  ans = data
  for(j in 1:ncol(data)){
    ans[,j] = col2cop(data[,j], denom_method)
  }
  return(ans)
}
