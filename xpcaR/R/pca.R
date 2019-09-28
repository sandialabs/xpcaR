#' @export
#' @title PCA
#' @param data matrix or data.frame
#' @param rank Rank of decomposition
#' @param tol Rank of tolerance used to assess convergence
#' @param rescale Should each column be scaled?
#' @param max_iter Maximum iterations
#' @param post_svd Should decomposition be post processed with svd?
#' @examples 
#' data = simProblem()$data
#' fit = pca(data, rank = 2)
pca <- function(data, rank, 
                tol = 0.1, 
                rescale = T, 
                max_iter = 500, 
                post_svd = T){
  data <- as.matrix(data)
  nRow = nrow(data)
  nCol = ncol(data)
  
  # Standardizing data
  col_means = colMeans(data, na.rm = TRUE)
  col_sds = rep(NA, nCol)
  this_sd = 1
  for(j in 1:nCol){
    if(rescale) this_sd <- samp_sd(data[,j])
    col_sds[j] = this_sd
    data[,j] <- (data[,j] - col_means[j]) / this_sd
  }
  
  # Matrix of logicals, indicated that data is not missing
  use = !is.na(data) 

  # Initializing A and B
  A = matrix(rnorm(nRow * rank), nrow = nRow)
  B = matrix(rnorm(nCol * rank), nrow = nCol)
  
  theta = A %*% t(B)
  loss = sum( (theta[use] - data[use])^2)
  err = tol + 1
  iter = 0
  
  while(tol < err & iter < max_iter){
    iter = iter + 1
    A = updateA(A, B, data, use)
    B = updateB(A, B, data, use)
    theta = A %*% t(B)
    loss_new = sum( (theta[use] - data[use])^2)
    err = loss - loss_new
    loss = loss_new
  }
  
  if(post_svd){
    svd_theta = svd(theta)
    A = svd_theta$u %*% diag(svd_theta$d)
    B = svd_theta$v
    A = A[,1:rank, drop = F]
    B = B[,1:rank, drop = F]
  }
  
  data_ests = theta
  for(j in 1:nCol){
    data_ests[,j] = theta[,j] * col_sds[j] + col_means[j]
  }
  
  rownames(A) = rownames(data)
  colnames(A) = paste0("PC", 1:rank)
  rownames(B) = colnames(data)
  colnames(B) = paste0("PC", 1:rank)
  dimnames(data_ests) = dimnames(data)
  
  ans = list(fittedEsts = data_ests, 
             theta = theta, 
             A = A, 
             B = B, 
             loss = loss, 
             iters = iter)
  return(ans)
}


# Computes the standard deviation of vector
# R's sd() function not used because that uses denominator n-1 instead of n
samp_sd = function(x){ sqrt( mean((x - mean(x, na.rm = T))^2, na.rm = T))}


# Simplified least squares for alternating least squares update
# R's lm spends way too much time data munging
subset_ls = function(X, y, use){
  y = subset(y, use)
  X = subset(X, use)
  nRow = nrow(X)
  nCol = ncol(X)
  if(nCol > nRow){
    ans = singular_ls(X, y)
    return(ans)
  }
  t_X = t(X)
  ans = as.numeric( solve(t_X %*% X, t_X %*% y) )
  return(ans)
}

# Due to missing data, sometimes the least squares solution
#  will be unidentified. This identifies the solution
# by fixing the unnecessary coefficients at 0. 
singular_ls = function(X, y){
  nRow = nrow(X)
  nCol = ncol(X)
  if(nRow == 0) return(rep(0, nCol))
  nDrop = nCol - nRow
  X_skinny = X[,1:(nCol - nDrop)]
  t_X_skinny = t(X_skinny)
  skinny_solution = as.numeric( solve(t_X_skinny %*% X_skinny,
                                      t_X_skinny %*% y) )
  ans = c(skinny_solution, rep(0, nDrop))
  return(ans)
}



# Updates all rows of A 
# "weights" are really just 1 if observed, 0 if missing
updateA = function(A, B, data, use){
  nRow = nrow(A)
  for(i in 1:nRow){
    this_y = data[i,]
    this_use = use[i,]
    new_row = subset_ls(B, this_y, this_use)
    A[i,] = new_row
  }
  return(A)
} 

# Updates all rows of B
updateB = function(A, B, data, use){
  nRow = nrow(B)
  for(i in 1:nRow){
    this_y = data[,i]
    this_use = use[,i]
    new_row = subset_ls(A, this_y, this_use)
    B[i,] = new_row
  }
  return(B)
} 
