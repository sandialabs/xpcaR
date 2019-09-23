# The mean in the case of exponential marginals 
# is not in closed form. Estimates using Monte Carlo estimation
MC_expMean <- function(MC = 1000, 
                       cop_mu, 
                       sigma, 
                       q_fun = qexp){
  means = NULL
  ses   = NULL
  if(length(sigma) == 1){ sigma = rep(sigma, length(cop_mu)) }
  for(i in seq_along(cop_mu)){
    cop_samps <- cop_mu[i] + rnorm(MC, sd = sigma[i])
    u_samps <- pnorm(cop_samps)
    x_samps <- q_fun(u_samps)
    means[i] = mean(x_samps)
    ses[i] = sd(x_samps) / sqrt(MC)
  }
  ans <- data.frame(means = means, se = ses)
  return(ans)
}

#' @export
simProblem <- function(nRow = 100, nCol = 100, 
                        sigma = 0.5, propBinary = 0, 
                        contError = "normal", rank = 3){
  if(sigma > 1 | sigma < 0) stop("sigma must be in [0,1]")
  # Step 1: Simulating V
  V_raw = matrix(rnorm(rank * nCol), nrow = nCol)
  VV_t_raw = V_raw %*% t(V_raw)
  VV_t = cov2cor(VV_t_raw) * (1 - sigma^2)
  
  # Step 2: Simulating Theta. Note that the random noise
  # has NOT been placed on top yet
  theta = mvtnorm::rmvnorm(nRow, sigma = VV_t)
  
  # Step 4: Setting binary columns
  isBinary = rep(F, nCol)
  isBinary[seq_len(round(nCol * propBinary))] = T
  
  # Step 5: Computing means for each column
  meanMatrix = matrix(nrow = nRow, ncol = nCol)
  for(j in 1:nCol){
    if(isBinary[j]){ 
      meanMatrix[,j] = pnorm(theta[,j], sd = sigma)
    }
    else{
      if(contError == "normal"){
        meanMatrix[,j] = theta[,j]
      }
      else if(contError == "exp"){
        thisColMean = MC_expMean(cop_mu = theta[,j], 
                                 sigma = sigma)
        meanMatrix[,j] = thisColMean$means
      }
    }
  }
  # Step 6: Simulating data
  copMatrix = theta + rnorm(nRow * nCol, sd = sigma)
  dataMatrix = matrix(nrow = nRow, ncol = nCol)
  for(j in 1:nCol){
    if(isBinary[j]){ dataMatrix[,j] = copMatrix[,j] > 0 }
    else{
      if(contError == "normal"){
        dataMatrix[,j] = copMatrix[,j]
      }
      else if(contError == "exp"){
        dataMatrix[,j] = qexp(pnorm(copMatrix[,j]))
      }
    }
  }
  ans = list(data = dataMatrix, 
             meanMat = meanMatrix)
  return(ans)
}

#' @export
simProblemOld = function(m = 100, 
                      n = 100, 
                      propBinary = .9, 
                      rank = 2, 
                      noiseVar = .5, 
                      contIsExp = F){

  # Simulating factors
  A = matrix(rnorm(m * rank), nrow = m)
  B = matrix(rnorm(n * rank), nrow = n)
  theta = A %*% t(B)
#  for(j in 1:n){ theta[,j] <- theta[,j] / sd(theta[,j])}
  meanMat = theta
  
  # Creating data 
  isBin = rep(FALSE, n)
  if(propBinary > 0){isBin[1:(n*propBinary)] <- TRUE}
  data = theta + sqrt(noiseVar) * rnorm(m*n)
  meanMat = theta
  
  # Making propBinary proportion of them binary
  data[,isBin] <- data[,isBin] > 0
  
  # Switching to exponential if specified
  if(contIsExp){
    shape = 5
    data[,!isBin] <- rgamma(length(theta[,!isBin]), 
                            shape = shape,
                            scale = exp(theta[,!isBin]) / shape )
    meanMat[,!isBin] <- exp(theta[,!isBin])
  }  
  meanMat[,isBin] <- pnorm(theta[,isBin]/noiseVar)

  ans = list(data = data, 
             meanMat = meanMat)
  return(ans)
}
