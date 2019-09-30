#' @title XPCA
#' @description Computes low-rank XPCA decomposition of data. 
#' Can be used with missing data. If all data is continuous, 
#' estimates will be very close to those provided by COCA. 
#' If marginal distribution of all columns is approximately 
#' Gaussian, both COCA + XPCA estimates will be very close 
#' to those provided by PCA. 
#' @param data Data matrix or data.frame
#' @param rank Rank of decomposition
#' @param tol Tolerance; if change in LLK is below this value, algorithm terminates
#' @param maxIters Maximum iterations
#' @param opt Optimization algorithm to use. Options are "LBFGS" or "AN" (Alternating Newton's)
#' @param reg_val Value used to regularize max intervals. Prevents Hauck-Donner effect. 
#' @param post_svd Should svd be used to standardize decomposition?
#' @param gridSize Number of points used to approximate mean function
#' @details Copula-based decomposition, 
#' where all data is treated as ordinal. 
#' 
#' For technical details, see:
#' 
#' C. Anderson-Bergman, T. G. Kolda, K. Kincher-Winoto. 
#' XPCA: Extending PCA for a 
#' Combination of Discrete and Continuous Variables. 
#' arXiv:1808.07510, 2018
#' 
#' @return 
#' A list with the following fields:
#' 
#'    - \code{A}: Low-rank representation of rows
#'    
#'    - \code{B}: Low-rank representation of columns
#'    
#'    - \code{fittedEsts}: Data estimates of low-rank 
#'    representation of full data matrix
#'    
#' Other fields are for internal use. 
#' @examples
#'data = simProblem()$data
#'fit = xpca(data, rank = 2)
#' @export
xpca = function(data, rank, tol = 0.1,
                maxIters = 1000, 
                opt = "AN", 
                reg_val = 0.5, 
                post_svd = T, 
                gridSize = nrow(data)/10){
  updateSigma = T
  data = as.matrix(data)
  # Construct [L,R] interval for each data point
  LR_list = raw2intervals(data, "m", reg_val)
  nRow = nrow(data)
  nCol = ncol(data)
  # Initial values for A, B matrix
  A = matrix(rnorm(rank*nRow, sd = 0.1), nrow = nRow)
  B = matrix(rnorm(rank*nCol, sd = 0.1), nrow = nCol)
  
  # Check if opt is one we recognize
  if (opt!="AN" &&  opt!="LBFGS") { 
      stop("opt argument not recognized. Options = LBFGS or AN") 
  }

  # Building XPCA optimizer object
  if(opt == "LBFGS"){
    xpca_opter = xpcaLBFGSClass(A, B, LR_list$L, LR_list$R, 
                           update_s = updateSigma)
    # Getting info required for L-BFGS-B algorithm
    startVals = xpca_opter$getVecVals()
    
    # Optimizing with L-BFGS-B
    ktot = length(startVals)
    lowerBnds <- rep(-4, ktot)
    upperBnds <- rep(4, ktot)
    lowerBnds[ktot] <- 0.01
    upperBnds[ktot] <- 1
     
    opt_res = nloptr::lbfgs(startVals, 
                      fn = xpca_opter$checkNewLossVec, 
                      gr = xpca_opter$flatDerv_fromVec,
                      lower = lowerBnds, upper = upperBnds,
                      control = list(ftol_abs = tol, 
                                     maxeval = maxIters)
                      )
  
    if(opt_res$convergence < 0){ 
      warning(paste("L-BFGS did not converge. Algorithm message:\n", opt_res$message), 
                    "\nTrying again with Alternating Newton's\n")
      opt = "AN"
        #return (NaN) # return NaN for timing purposes
    }
    else {
        # Extracting vector optimal values 
        flat_vals = opt_res$par
        # Plugging back into xpca_opter to reshape
        xpca_opter$setVecVals(flat_vals)
        final_loss = opt_res$value
        iters_used = opt_res$iter
    }
  }
  if(opt == "AN"){
    xpca_opter = xpcaAltNewtClass(A, B, 
                                  LR_list$L, LR_list$R, 
                                  update_s = updateSigma)
    opt_res <- xpca_opter$alt_newtons(tol = tol, 
                           max_iters = maxIters, 
                           verbose = F)
    
    final_loss = opt_res["loss"]
    iters_used = opt_res["iters"]
  }
  
  # Extracting components
  A = xpca_opter$A
  B = xpca_opter$B
  theta = A %*% t(B)
  if(post_svd){
    svd_theta = svd(theta)
    A = svd_theta$u %*% diag(svd_theta$d)
    B = svd_theta$v
    A = A[,1:rank, drop = F]
    B = B[,1:rank, drop = F]
  }
  sigma = xpca_opter$sigma
  col_sds = NULL
  if(!updateSigma){
    for(j in 1:ncol(data)){ col_sds[j] = xpca_get_col_sd(j, theta, LR_list$L, LR_list$R) }
   sigma = mean(col_sds)
  }
  fittedMeans = matrix(nrow = nrow(data), 
                       ncol = ncol(data))
  cdfs = makeAllEDFs(data)
  for(j in 1:ncol(data)){
    fittedMeans[,j] = computeColMean_approx(j, sigma = sigma, 
                                            z_matrix = theta, 
                                            cdfs = cdfs, 
                                            gridSize = gridSize)
  }
  rownames(A) = rownames(data)
  colnames(A) = paste0("PC", 1:rank)
  rownames(B) = colnames(data)
  colnames(B) = paste0("PC", 1:rank)
  colnames(fittedMeans) = colnames(data)
  ans =  list(A = A, B = B, 
              theta = theta, 
              sigma = sigma,
              fittedEsts = fittedMeans, 
              loss = final_loss, 
              iters = iters_used, 
              columnCDFs = cdfs)
  class(ans) <- "xpca_class"
  return(ans)
}


### PREPROCESSING STEP

# Takes in one column of raw data 
# Outputs the left and right interval to be passed to optimization step 
col2intervals = function(vals, denom_method = "m", reg_val){
  m = sum(!is.na(vals))
  
  # "m" method is used for XPCA, "m+1" is used for COCA
  # Slightly different ways of calculating Empirical Distribution Function
  if(denom_method == "m"){denom = m}
  else if(denom_method == "m+1"){denom = m+1}
  else{stop("denom_method invalid. Valid options: 'm' or 'm+1'")}
  L = (rank(vals, ties.method = "min", na.last = "keep")-1)/denom
  R = rank(vals, ties.method = "max", na.last = "keep")/denom
  
  # As a form of regularization, will pull in the outmost intervals
  if(denom_method == "m" ){
    minVal = min(vals, na.rm = T)
    maxVal = max(vals, na.rm = T)
    smallestValInds = which(vals == minVal)
    largestValInds = which(vals == maxVal)
    L[smallestValInds] = reg_val/denom
    R[largestValInds] = 1 - reg_val/denom
  }
  ans = cbind(L, R)
  return(ans)
}

# Takes in full data matrix
# Outputs full L, R matrices 
raw2intervals = function(data, denom_method = "m", reg_val){
  data = as.matrix(data)
  # Creating Left and Right side of intervals for data
  L = data * 0 - 1
  R = data * 0 - 1
  for(j in 1:ncol(data)){
    intervals = col2intervals(data[,j], denom_method, reg_val)
    L[,j] = intervals[,1]
    R[,j] = intervals[,2]
  }
  # L[is.na(L)] = 0
  # R[is.na(R)] = 1
  
  L[is.na(L)] = min(L, na.rm = T)
  R[is.na(R)] = max(L, na.rm = T)
  
  ans = list(L, R)
  names(ans) = c("L", "R")
  return(ans)
}


### OPTIMIZATION STEP
# Exported to xpcaAltNewtClass or xpcaLBFGSClass


### POST PROCESSING STEP

getFrequencies = function(vals){
  # table(x) gives you the count of each unique value in x
  counts = table(vals)
  # names(table(x)) is the actual values themselves. Need to convert back to numbers
  vals = as.numeric(names(counts))
  names(counts) = NULL
  ans = cbind(vals, counts)
}

# Takes in one column of data, outputs the empirical distribution function
makeEDF <- function(vals, denom_method = "m"){
  m <- sum(!is.na(vals) ) 
  if(denom_method == "m") denom = m
  else if(denom_method == "m+1") denom = m+1
  else stop("denom_method should be 'm' or 'm+1'")
  obs_vals = vals[!is.na(vals)]
  freq_counts = getFrequencies(obs_vals)
  cdf_x = freq_counts[,1] 
  cdf_p = cumsum(freq_counts[,2]) / denom   
  ans = list(x = cdf_x, p = cdf_p)
  return(ans)
}

# Takes in full data, returns list of CDF's
makeAllEDFs = function(data, denom_method = "m"){
  ans = list()
  for(j in 1:ncol(data)){
    ans[[j]] = makeEDF(data[,j], denom_method)
  }
  return(ans)
}

### COMPUTING ESTIMATED SIGMAS IN POSTPROCESSING

# Compute loglikelihood for sigma based point estimate of mean
# To be passed to generic optimization method 
# sigma is input to be optimized, mean_vec is a column of means (i.e. theta[,j])
# L_vec and R_vec are column vectors of the L, R matrix 
icen_sd_nllk <- function(sigma, mean_vec, L_vec, R_vec){
  ans <- 0
  ans <- log(sigma) 
  icenProbs <- pnorm(q = R_vec, 
                     mean = mean_vec, 
                     sd = sigma) - 
               pnorm(q = L_vec, 
                     mean = mean_vec, 
                     sd = sigma)
    ans <- ans + sum( log(icenProbs) )
  # Note: R's optimization methods fail when non-finite value is returned
  # Below is a workaround. 
  if(!is.finite(ans)) ans <- -10^15
  return(-ans)
}

# Compute conditional standard deviation conditional on theta estimates
icen_sd_mle <- function(mean_vec, L_vec, R_vec, rng = c(0,1)){
  # Need to transform L, R, which is in [0,1], to z-scores
  L_z = qnorm(L_vec)
  R_z = qnorm(R_vec)
  optRes <- optimize(icen_sd_nllk, interval = rng,
                     mean_vec = mean_vec, 
                     L_vec = L_z, 
                     R_vec = R_z)
  return(optRes$minimum)
}

# extracts the estimated conditional column sd's
xpca_get_col_sd <- function(col_index, z_est, 
                       leftMatrix, rightMatrix){
  this_z     <- z_est[,col_index]
  this_left  <- leftMatrix[, col_index]
  this_right <- rightMatrix[,col_index]

  ans <- icen_sd_mle(this_z, this_left, this_right)
  return(ans)
}

### COMPUTING MEAN ESTIMATES IN POST PROCESSING


conditionalProbs <- function(this_z, this_s, cdf_z){
  # Recentering and rescaling conditional on latent mean + sigma
  recentered_z = (cdf_z - this_z) / this_s
  # Computing conditional cdf
  cdf_p = pnorm(recentered_z)
  # Computing probabilities of each entry
  k = length(cdf_p)
  if(cdf_p[k] != 1){ cdf_p = cdf_p / cdf_p[k] }
  ans = cdf_p - c(0, cdf_p[-k])
  return(ans)
}

# Computes the mean for a given z value, estimated standard deviation and cdf
computeMean <- function(this_z, this_s, cdf_x, cdf_z){
  probs = conditionalProbs(this_z, this_s, cdf_z)
  ans          <- sum(probs * cdf_x)
  return(ans)
}



# Computes approximated means by computing mean across a grid 
# of z values and then using linear interpolation (i.e. approxfun)
computeColMean_approx <- function(col_index, sigma, 
                                  z_matrix, cdfs, 
                                  approxMaker = approxfun, 
                                  gridSize = 100){
  # Note: approxfun is an R function that creates a linear interpolator for a set of values

  # Extracting necessary information for this column
  these_zs = z_matrix[,col_index]
  this_cdf = cdfs[[col_index]]
  
  # First computing means over a grid of inputs, ranging from min to max z-value 
  sample_zs = quantile(these_zs, probs =  0:gridSize / gridSize)
  sampleMeans = NULL
  for(i in seq_along(sample_zs)){
    sampleMeans[i] = computeMean(sample_zs[i], 
                                  this_s = sigma, 
                                  cdf_x = this_cdf$x, 
                                  cdf_z = qnorm(this_cdf$p))
  }
  isNA <- is.na(sampleMeans)
  if(sum(isNA) > 0){
    warning(paste("Column", col_index, 'fitted NAs during approximation step'))
    sampleMeans <- sampleMeans[!isNA]
    sample_zs <- sample_zs[!isNA]
  }
  if(length(sampleMeans) == 0){
    ans <- rep(NA, length(these_zs))
    return(ans)
  }
  if(all(sampleMeans == sampleMeans[1])){
    warning(paste("Note: Column", col_index, 'is fitted as constant'))
    ans <- rep(sampleMeans[1], length(these_zs))
    ans <- as.numeric(ans)
    return(ans)
  }
  fitter <- approxMaker(x = as.numeric(sample_zs), 
                        y = as.numeric(sampleMeans), 
                        rule = 2)
  theseMeans <- fitter(these_zs)
  return(theseMeans)
}

## USER UTILITY FUNCTIONS
#  Currently not exported

getConditionalPMF <- function(i, j, xpca_fit){
  if(!is(xpca_fit, "xpca_class")) stop("xpca_fit must be of class xpca_class")
  s = xpca_fit$s
  cdfs <- xpca_fit$columnCDFs
  theta <- xpca_fit$theta
  this_z = theta[i,j]
  this_cdf <- cdfs[[j]]
  probs = conditionalProbs(this_z, s,
                           cdf_z = qnorm(this_cdf$p))
  ans = cbind(this_cdf$x, probs)
  colnames(ans) = c("x", "p")
  rownames(ans) <- NULL
  return(ans)
}

getPredictionInterval <- function(i, j, 
                                  xpca_fit, 
                                  predLevel = .95){
  interval_probs = (1 - predLevel)/2 * c(1, -1) + c(0, 1)
  pmf = getConditionalPMF(i, j, xpca_fit)
  cum_p = cumsum(pmf$p)
  if(cum_p[1] > interval_probs[1]) lower_ind = 1
  else{ lower_ind = max(which(cum_p < interval_probs[1])) }
  
  upper_ind = min(which(cum_p >= interval_probs[2]))
  ans = c(pmf$x[lower_ind], 
          pmf$x[upper_ind])
  names(ans) = c(paste0("lower_", interval_probs[1] ),
                 paste0("upper_", interval_probs[2]))
                 
  return(ans)
}
