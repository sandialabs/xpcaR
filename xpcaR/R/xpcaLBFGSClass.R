# xpcaLBFGSClass is a class for the optimization step ONLY
xpcaLBFGSClass = setRefClass("xpcaLBFGSClass", 
                        fields = c("L", "R", 
                                   "A", "B", 
                                   "theta", "sigma", 
                                   "update_sigma"), 
                        methods = c("initialize",   
                                    # constructor
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
                                    "getVecVals"
                                    # Takes current values of A,B and sigma and flattens
                        ))

## SETUP, GETTERS AND SETTERS
xpcaLBFGSClass$methods(
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
    update_sigma <<- update_s
  }
)



xpcaLBFGSClass$methods(
  setVecVals = function(all_vec){
    A_length = length(A)
    B_length = length(B)
    if(length(all_vec) != A_length + B_length + 1) 
      stop("input to checkNewLossVec wrong length")
    A_vals = all_vec[1:A_length]
    B_vals = all_vec[A_length + 1:B_length]
    sigma_val = all_vec[A_length + B_length + 1]
    A_new = matrix(A_vals, nrow = nrow(A))
    B_new = matrix(B_vals, nrow = nrow(B))
    A <<- A_new
    B <<- B_new
    sigma <<- sigma_val
  }
)

xpcaLBFGSClass$methods(
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

xpcaLBFGSClass$methods(
  checkNewLoss = function(A_new, B_new, sigma_new){
      # Negative Log Likelihood for XPCA
    theta_new = A_new %*% t(B_new)
    L_adj = (L - theta_new)/sigma_new
    R_adj = (R - theta_new)/sigma_new
    probs = pnorm(R_adj) - pnorm(L_adj)
    ans = -sum(log(probs))
    return(ans)
  }
)

xpcaLBFGSClass$methods(
  checkNewLossVec = function(all_vec){
    A_length = length(A)
    B_length = length(B)
    if(length(all_vec) != A_length + B_length + 1) 
      stop("input to checkNewLossVec wrong length")
    A_vals = all_vec[1:A_length]
    B_vals = all_vec[A_length + 1:B_length]
    sigma_val = all_vec[A_length + B_length + 1]
    A_new = matrix(A_vals, nrow = nrow(A))
    B_new = matrix(B_vals, nrow = nrow(B))
    ans = checkNewLoss(A_new, B_new, sigma_val)
    return(ans)
  }
)

xpcaLBFGSClass$methods(
  flatDerv_fromVec = function(input){
    A_length = length(A)
    B_length = length(B)
    if(length(input) != A_length + B_length + 1) 
      stop("inpute length incorrect: internal error")
    A_vals = input[1:A_length]
    B_vals = input[A_length + 1:B_length]
    sigma_new = input[A_length + B_length + 1]
    A_new = matrix(A_vals, nrow = nrow(A))
    B_new = matrix(B_vals, nrow = nrow(B))
    ans = flatDerv(A_new, B_new, sigma_new)
    return(ans)
  }
)
xpcaLBFGSClass$methods(
  flatDerv = function(A_new, B_new, sigma_new){
    theta_new = A_new %*% t(B_new)
    L_adj = (L - theta_new)/sigma_new
    R_adj = (R - theta_new)/sigma_new
    probs = pnorm(R_adj) - pnorm(L_adj)
    phi_l = dnorm(L_adj) 
    phi_r = dnorm(R_adj) 
    # If L = -Inf, we need L_adj * phi(L_adj) == 0
    # Same for R_adj == Inf
    R_adj[R_adj == Inf] = 0
    L_adj[L_adj == -Inf] = 0

    d1_base = -(phi_l - phi_r)/( probs * sigma_new)
    A_dervs = d1_base %*% B_new
    B_dervs = t(d1_base) %*% A_new
    if(update_sigma){ 
      sigma_derv = sum( (phi_r * R_adj - phi_l * L_adj) / probs ) / sigma_new
    }
    else{ sigma_derv = 0 }
    ans = c(as.numeric(A_dervs),
            as.numeric(B_dervs), 
            sigma_derv)
    return(ans)
  }
)
