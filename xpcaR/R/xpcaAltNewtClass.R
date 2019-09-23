# xpcaAltNewtClass is a class for the optimization step ONLY
xpcaAltNewtClass = setRefClass("xpcaAltNewtClass", 
                        fields = c("A", "B", 
                                   "L", "R", 
                                   "theta", 
                                   "sigma", 
                                   "currentLoss", 
                                   "updateSigma"), 
                        methods = c("initialize",   
                                    # constructor
                                    "computeProbs", 
                                    # computes probability for a vector of left side of interval, 
                                    # right side of interval and corresponding theta's
                                    "computeBaseDervs", 
                                    # computes derivatives *with respect to theta*
                                    "flatDerv",
                                    # computes a vector with derivatives of A, B and sigma 
                                    "flatDerv_fromVec",
                                    # same as above, but takes in flat vector
                                    "checkNewLoss",
                                    # computes, without setting parameters, new log-likelihood 
                                    "checkNewLossVec",
                                    # Takes in a flat vector containing values of A, B and sigma
                                    # and checks Loss without setting values
                                    "setVecVals",
                                    # Takes in flat vector and sets A, B, and sigma
                                    "getVecVals",
                                    # Takes current values of A,B and sigma and flattens
                                    "newtonShift",
                                    # computes the newton proposed change in parameters
                                    "halfStep",
                                    # Uses half-stepping to ensure decrease in loss when updating A & B
                                    "updateA",
                                    # updates all the rows of A with newton's method
                                    "updateB",
                                    # updates all the rows of B with newton's method
                                    "alt_newtons" 
                                    # runs alternating Newton's method algorithm
                        ))

## SETUP, GETTERS AND SETTERS
xpcaAltNewtClass$methods(
  initialize = function(A, B, L, R, update_s){
    # R note: inside a reference class method, "<<-" means set field
    # If you wrote "myField <- A", this would make a local object called "myField"
    # rather than copy into the field called myField
    A <<- A
    B <<- B
    theta <<- A %*% t(B)
    # L and R should probabilities, which are transformed to z-scores
    L <<- qnorm(L) 
    R <<- qnorm(R) 
    sigma <<- 1
    updateSigma <<- update_s
  }
)

expit = function(x) exp(x) / (1 + exp(x))

xpcaAltNewtClass$methods(
  setVecVals = function(all_vec){
    A_length = length(A)
    B_length = length(B)
    if(length(all_vec) != A_length + B_length + 1) 
      stop("input to checkNewLossVec wrong length")
    A_vals = all_vec[1:A_length]
    B_vals = all_vec[A_length + 1:B_length]
    lgt_sigma_val = all_vec[A_length + B_length + 1]
    A_new = matrix(A_vals, nrow = nrow(A))
    B_new = matrix(B_vals, nrow = nrow(B))
    A <<- A_new
    B <<- B_new
    sigma <<- sigma_val
  }
)

xpcaAltNewtClass$methods(
  getVecVals = function(){
    ans = c(as.numeric(A), 
            as.numeric(B), 
            sigma)
    return(ans)
  }
)

## L-BFGS SPECIFIC FUNCTIONS
#  These functions should take in 
# vectors of parameters rather than cached values

xpcaAltNewtClass$methods(
  checkNewLoss = function(A_new, B_new, sigma_new){
#    sigma_new = expit(lgt_sigma_new)
    theta_new = A_new %*% t(B_new)
    L_adj = (L - theta_new)/sigma_new
    R_adj = (R - theta_new)/sigma_new
    probs = pnorm(R_adj) - pnorm(L_adj)
    return(-sum(log(probs)))
  }
)

xpcaAltNewtClass$methods(
  checkNewLossVec = function(all_vec){
    A_length = length(A)
    B_length = length(B)
    if(length(all_vec) != A_length + B_length + 1) 
      stop("input to checkNewLossVec wrong length")
    A_vals = all_vec[1:A_length]
    B_vals = all_vec[A_length + 1:B_length]
    lgt_sigma_val = all_vec[A_length + B_length + 1]
    A_new = matrix(A_vals, nrow = nrow(A))
    B_new = matrix(B_vals, nrow = nrow(B))
    ans = checkNewLoss(A_new, B_new, lgt_sigma_val)
    return(ans)
  }
)

xpcaAltNewtClass$methods(
  flatDerv_fromVec = function(input){
    A_length = length(A)
    B_length = length(B)
    if(length(input) != A_length + B_length + 1) 
      stop("inpute length incorrect: internal error")
    A_vals = input[1:A_length]
    B_vals = input[A_length + 1:B_length]
    lgt_sigma_new = input[A_length + B_length + 1]
    A_new = matrix(A_vals, nrow = nrow(A))
    B_new = matrix(B_vals, nrow = nrow(B))
    ans = flatDerv(A_new, B_new, lgt_sigma_new)
    return(ans)
  }
)
xpcaAltNewtClass$methods(
  flatDerv = function(A_new, B_new, sigma_new){
    theta_new = A_new %*% t(B_new)
    L_adj = (L - theta_new)/sigma_new
    R_adj = (R - theta_new)/sigma_new
    probs = pnorm(R_adj) - pnorm(L_adj)
    phi_l = dnorm(L_adj) 
    phi_r = dnorm(R_adj) 
    d1_base = -(phi_l - phi_r)/( probs * sigma_new)
    A_dervs = d1_base %*% B_new
    B_dervs = t(d1_base) %*% A_new
    sigma_derv = sum( (phi_r * R_adj - phi_l * L_adj) / probs ) / sigma_new
    ans = c(as.numeric(A_dervs),
            as.numeric(B_dervs), 
            sigma_derv)
    return(ans)
  }
)

## METHODS FOR ALTERNATING OPTIMIZATION

# Base function for computing probabilities
xpcaAltNewtClass$methods(
  computeProbs = function(L_vec, R_vec, theta_vec){
    # computes probability for a vector of left and right side values, 
    # with corresponding theta's
    R_adj = (R_vec - theta_vec) / sigma
    L_adj = (L_vec - theta_vec) / sigma
    probs = pnorm(R_adj) - pnorm(L_adj)
    return(probs)
  })


xpcaAltNewtClass$methods(
  full_loss = function(){
    # computes loss across all values
    # Note that this is the *negative* log likelihood (convex)
    probs = computeProbs(L, R, theta)
    ans = -sum(log(probs))
    if(!is.finite(ans)) ans = Inf
    return(ans)
  })

xpcaAltNewtClass$methods(
  computeBaseDervs = function(L_vec, R_vec, theta_vec){
    # Computes derivatives for a vector of left, right side values and theta
    # Computes derivatives with respect to *theta*
    L_adj = (L_vec - theta_vec)/sigma
    R_adj = (R_vec - theta_vec)/sigma
    probs = pnorm(R_adj) - pnorm(L_adj)
    phi_l = dnorm(L_adj)/sigma
    phi_r = dnorm(R_adj)/sigma
    phi_l_prime = -phi_l * L_adj
    phi_r_prime = -phi_r * R_adj
    
    phi_l_prime[L_vec == -Inf] = 0
    phi_r_prime[R_vec == Inf] = 0
    
    # d1 = first derivative
    # d2 = second derivative
    d1 = (phi_l - phi_r) / probs
    d2 = (phi_l_prime - phi_r_prime) / probs + d1 * d1
    
    ans = cbind(d1, d2)
    return(ans)
  })


xpcaAltNewtClass$methods(
  # computes proposed shift from Newton's method 
  # Inpute d1_vec, d2_vec is first and second derivatives *with respect to theta*
  # X is either the A or B matrix (if we are updating a row A, this is all of B and vice versa)
  newtonShift = function(d1_vec, d2_vec, X){
    # Computing derivatives with respect to A or B matrix
    if(all(d2_vec > 0) ) 
      X_adj = X * sqrt(d2_vec)
    else
      X_adj = X 
    hess = t(X_adj) %*% X_adj
    d1 = t(X) %*% d1_vec
    # Solving for the proposal step  
    ans = solve(hess, d1)
    if(any(!is.finite(ans))) return(rep(0, length(ans)))
    return(ans)
  })

xpcaAltNewtClass$methods(
  halfStep = function(L_vec, R_vec, theta_vec, beta, X, propDelta){
    # This method is used to ensure that the loss decreases at each step
    # Computing starting loss from this row
    startingLoss = -sum(log(computeProbs(L_vec, R_vec, theta_vec)))
    # Computing new loss at proposal 
    new_beta = beta + propDelta
    new_theta = X %*% new_beta
    newLoss = -sum(log(computeProbs(L_vec, R_vec, new_theta)))
    if(!is.finite(newLoss)) newLoss <- Inf
    # Checking to see if loss improved
    ok = newLoss < startingLoss
    # If not, we will try several half steps
    tries = 0
    while( (!ok) & (tries < 5) ){
      tries = tries + 1
      propDelta = propDelta / 2 
      new_beta = beta + propDelta
      new_theta = X %*% new_beta
      newLoss = -sum(log(computeProbs(L_vec, R_vec, new_theta)))
      if(!is.finite(newLoss)) newLoss <- Inf
      ok = newLoss <= startingLoss
    }
    # If we found a good solution, ok will be true and we will update the row
    # If not, we will not update the row
    if(ok){ 
      currentLoss <<- currentLoss + (newLoss - startingLoss)
      return(new_beta) 
    }
    else { return(beta) }
  }
)



xpcaAltNewtClass$methods(
  updateA = function(){
    # Updates all of A matrix
    nRows = nrow(A)
    # Placeholder for A. I *think* if we updated a row of the field A, 
    # it would trigger a deep copy of all of A? Not actually sure about that...
    temp_A = A
    for(i in 1:nRows){
      # Grabbing the necessary info from this row
      this_L = L[i,]
      this_R = R[i,]
      this_theta = theta[i,]
      # Computing the derivatives with respect to theta
      baseDervs = computeBaseDervs(this_L, this_R, this_theta)
      # Computing the proposed shift from Newton's method
      this_delta = newtonShift(d1_vec = baseDervs[,1], d2_vec = baseDervs[,2], B)
      # Double checking that it does decrease loss
      newRowVals = halfStep(this_L, this_R, this_theta, A[i,], B, this_delta)
      # Updating row i
      temp_A[i,] = newRowVals
    }
    # Updating the field A
    A <<- temp_A
    theta <<- A %*% t(B)
  }
)

xpcaAltNewtClass$methods(
  updateB = function(){
    # Updates all of B matrix; see updateA for more thorough commenting
    # Also note that while updateB updates rows of B,
    # it needs columns of L,R and theta, unlike updateA (needs rows)
    nCols = nrow(B)
    temp_B = B
    for(j in 1:nCols){
      this_L = L[,j]
      this_R = R[,j]
      this_theta = theta[,j]
      baseDervs = computeBaseDervs(this_L, this_R, this_theta)
      this_delta = newtonShift(d1_vec = baseDervs[,1], d2_vec = baseDervs[,2], A)
      newRowVals = halfStep(this_L, this_R, this_theta, B[j,], A, this_delta)
      # Updating row j
      temp_B[j,] = newRowVals
    }
    B <<- temp_B
    theta <<- A %*% t(B)
  }
)

xpcaAltNewtClass$methods(
  update_sigma = function(verbose){
    # Performs Newton's method on sigma
    startingLoss = currentLoss
#    theta <<- A %*% t(B)
    # Computing basic parts that go into derivatives
    L_adj = (L - theta)/sigma
    R_adj = (R - theta)/sigma
    phi_l = dnorm(L_adj)
    phi_r = dnorm(R_adj)
    probs = pnorm(R_adj) - pnorm(L_adj)
    
    # The limit of R_adj * dnorm(R_adj) -> 0 as R_adj -> Inf
    # Same for L_Adj
    # Need to replace these values so that R_adj * dnorm(R_adj) = 0
    L_adj[L_adj == -Inf] <- 0
    R_adj[R_adj ==  Inf] <- 0
    # Computing 1st and 2nd derivatives
    phi_L_adj_diff  = phi_l * L_adj - phi_r * R_adj
    d1 = -1 / sigma * sum( (phi_L_adj_diff) / probs )
    d2_part1 = (phi_L_adj_diff)^2
    d2_part2 = (R_adj^3 * phi_r - L_adj ^3 * phi_l)
    d2_part3 =  2 * (-phi_L_adj_diff)
    d2 = sigma^(-2) * sum(probs^(-2) * 
                            (d2_part1 + probs * (d2_part2 + d2_part3)))

    # Computing proposal step
    # If function is locally convex, use Newton
    if(d2 > 0) delta = -d1 / d2
    else delta = -sign(d1) / 10
    # We want to force sigma to be in [0.001, 1]
    # So we force steps size to respect bounds
    delta = max(-sigma, delta)
    delta = min(1 - sigma, delta)
    
    # Updating sigma and saving old value for line search
    old_sigma = sigma
    sigma <<- old_sigma + delta
    # Checking that loss decreased and reducing step size if it does not
    currentLoss <<- full_loss()
    steps = 0
    while(currentLoss > startingLoss & steps < 4){
      steps <- steps + 1
      delta <- delta / 4
      sigma <<- old_sigma + delta 
      currentLoss <<- full_loss()
    }
    # If line search failed, resetting the sigma
    if(currentLoss > startingLoss){ 
      if(verbose){
        cat("Note: Newton's method for sigma failed.",
                                  "Delta = ", round( delta, 5), 
                                  "d1 = ", d1, "d2 = ", d2, "\n")
      }
      sigma <<- old_sigma
      currentLoss <<- startingLoss
    }
  }
)

xpcaAltNewtClass$methods(
  alt_newtons = function(tol = 0.1, max_iters = 100, 
                         verbose = T, updateSigma = T){
    # Alternating Newton's method
    err = tol + 1
    currentLoss <<- full_loss()
    previousLoss = currentLoss
    iter = 0
    while(iter < max_iters & err > tol){
      iter = iter + 1
      if(verbose) cat("Loss = ", currentLoss, "\n")
      updateA()
      updateB()
      if(updateSigma){ update_sigma(verbose = verbose) }
      err = previousLoss - currentLoss
      previousLoss = currentLoss
    }
    if(verbose) cat("Loss = ", currentLoss, "\n", "Total iterations = ", iter, "\n")
    ans = c(currentLoss, iter)
    names(ans) = c("loss", "iters")
    return(ans)
  }
)
