sd_mle = function(x, na.rm = TRUE){
  x_bar = mean(x, na.rm = na.rm)
  x2_bar = mean(x^2, na.rm = na.rm)
  ans = sqrt(x2_bar - x_bar^2)
  return(ans)
}

#' @title Rescaled Mean Squared Error for xpca_cv
#' @param x_true Column of data
#' @param x_est Column of estimates
#' @param col_index Current column index
#' 
#' @description Example function for \code{xpca_cv}. 
#' Column by column, the original data and the 
#' estimates **only for the data that was left out** will be
#' fed into this function. As such, this function must be able to ignore
#' missing data
#' @export
cv_rmse = function(x_true, x_est, col_index = 1){
  ans <- mean( (x_true - x_est)^2, na.rm = T)
  ans <- ans / sd_mle(x_true, na.rm = T)^2
  return(ans)
}

#' @title Cross validation for MXD methods
#' @param data Matrix or data frame
#' @param folds Number of cross-validation folds
#' @param methods Decomposition methods used 
#' @param rank Rank of decomposition
#' @param lossFun Function for computing column by column loss
#' @examples 
#' data = simProblem()$data
#' xpca_cv(data, rank = 2)
#' @export
xpca_cv <- function(data, 
                    folds = 4, 
                    methods = c("pca", "xpca", "coca", "col_mean"), 
                    rank  = 1, 
                    lossFun = cv_rmse,
                    per.column=NA){
  
  cnames = colnames(data)
  data_mat = as.matrix(data)

  # Determining which data is observed
  # Get obs
  if (!is.na(per.column)) {
      m = nrow(data)
      # Get the rows in which that column has no missing
      is_obs = which(!is.na(data_mat[,per.column])) 
      # Convert to indices
      is_obs = is_obs + ((per.column-1)*m)
  }
  else {
      is_obs = which(!is.na(data_mat)) 
  }
  n_obs = length(is_obs)

  # Setting up random partitions
  n_perPartition = n_obs / folds
  cur_obs = is_obs
  partitions <- list()
  for(i in 1:(folds-1)){
    these_indices = sample(1:length(cur_obs), n_perPartition)
    partitions[[i]] <- cur_obs[these_indices]
    cur_obs <- cur_obs[-these_indices]
  }
  partitions[[folds]] <- cur_obs
  
  # Creating matrices that will be populated with estimates. 
  xpca_ests = data_mat + NA
  coca_ests = data_mat + NA
  pca_ests  = data_mat + NA
  col_mean_ests = data_mat + NA
  
  # Partitioning data into folds, computing decompositions 
  # and filling in ests matrices
  for(i in 1:folds){
    this_data = data_mat
    these_indices = partitions[[i]]
    this_data[these_indices] = NA
    if(any(methods == "xpca")){
      fit = xpca(this_data, rank = rank)
      xpca_ests[these_indices] = fit$meanEsts[these_indices]
    }
    if(any(methods == "pca")){
      fit = pca(this_data, rank = rank)
      pca_ests[these_indices] = fit$meanEsts[these_indices]
    }
    if(any(methods == "coca")){
      fit = coca(this_data, rank = rank)
      coca_ests[these_indices] = fit$meanEsts[these_indices]
    }
    if(any(methods == "col_mean")){
      these_col_mean_vec = colMeans(this_data, na.rm = TRUE)
      these_col_mean_mat = matrix(these_col_mean_vec, 
                                  byrow = TRUE, 
                                  nrow = nrow(data), 
                                  ncol = ncol(data))
      col_mean_ests[these_indices] = these_col_mean_mat[these_indices]
    }
  }
  # Creating summary stat for each method and saving full ests
  summary_vec = NULL
  summary_vec_names = NULL
  full_ests = list()
  colLoss_list = list()
  
  if(any(methods == "xpca")){
    loss_vec = NULL
    for(j in 1:ncol(data_mat)){ loss_vec[j] = lossFun( data_mat[,j], xpca_ests[,j], j)}
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec) )
    summary_vec_names = c(summary_vec_names, "xpca")
    colLoss_list[["xpca"]] = loss_vec
    full_ests[["xpca"]] = xpca_ests
  }
  if(any(methods == "pca")){
    loss_vec = NULL
    for(j in 1:ncol(data_mat)){ loss_vec[j] =  lossFun(data_mat[,j], pca_ests[,j], j)}
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "pca")
    colLoss_list[["pca"]] = loss_vec
    full_ests[["pca"]] = pca_ests
  }
  if(any(methods == "coca")){
    loss_vec = NULL
    for(j in 1:ncol(data_mat)){ loss_vec[j] = lossFun(data_mat[,j], coca_ests[,j], j) }
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "coca")
    colLoss_list[["coca"]] = loss_vec
    full_ests[["coca"]] = coca_ests
  }
  if(any(methods == "col_mean")){
    loss_vec = NULL
    for(j in 1:ncol(data_mat)){ loss_vec[j] = lossFun(data_mat[,j], col_mean_ests[,j], j) }
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "col_mean")
    colLoss_list[["col_mean"]] = loss_vec
    full_ests[["col_mean"]] = col_mean_ests
  }
  # Prepping output
  colLossMat = NULL
  for(i in seq_along(colLoss_list)){
    colLossMat <- cbind(colLossMat, colLoss_list[[i]])
  }
  colnames(colLossMat) = names(colLoss_list)
  names(summary_vec) = summary_vec_names
  ans = list(error_stats = summary_vec, 
             colLosses = as.data.frame(colLossMat),
             fullEsts = full_ests)
  return(ans)
}

#' @title In Sample Error for MXD methods
#' @param data Matrix or data frame
#' @param methods Decomposition methods used 
#' @param rank Rank of decomposition
#' @param lossFun Function for computing column by column loss
#' @examples 
#' data = simProblem()$data
#' xpca_ise(data, rank = 2)
#' @export
xpca_ise = function(data, 
                   methods = c("pca", "xpca", "coca", "col_mean"), 
                   rank  = 1, 
                   lossFun = cv_rmse){
  summary_vec = NULL
  summary_vec_names = NULL
  colLoss_list = list()
  fullEsts = list()
  cnames <- colnames(data)
  if(any(methods == "xpca")){
    fit <- xpca(data, rank)
    ests <- fit$meanEsts
    loss_vec <- NULL
    for(j in 1:ncol(data)){ loss_vec[j] = lossFun(data[,j], ests[,j], j)}
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "xpca")
    colLoss_list[["xpca"]] = loss_vec
    fullEsts[["xpca"]] = ests
  }
  
  if(any(methods == "pca")){
    fit <- pca(data, rank)
    ests <- fit$meanEsts
    loss_vec <- NULL
    for(j in 1:ncol(data)){ loss_vec[j] = lossFun(data[,j], ests[,j], j)}
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "pca")
    colLoss_list[["pca"]] = loss_vec
    fullEsts[["pca"]] = ests
  }  
  
  if(any(methods == "coca")){
    fit <- coca(data, rank)
    ests <- fit$meanEsts
    loss_vec <- NULL
    for(j in 1:ncol(data)){ loss_vec[j] = lossFun(data[,j], ests[,j], j)}
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "coca")
    colLoss_list[["coca"]] = loss_vec
    fullEsts[["coca"]] = ests
  }
  
  if(any(methods == "col_mean")){
    ests_vec <- colMeans(data, na.rm = TRUE)
    ests = matrix(ests_vec, 
                      byrow = TRUE, 
                      nrow = nrow(data), 
                      ncol = ncol(data))
    loss_vec <- NULL
    for(j in 1:ncol(data)){ loss_vec[j] = lossFun(data[,j], ests[,j], j)}
    names(loss_vec) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec))
    summary_vec_names = c(summary_vec_names, "col_mean")
    colLoss_list[["col_mean"]] = loss_vec
    fullEsts[["col_mean"]] = ests
  }  
  
  names(summary_vec) = summary_vec_names
  ans = list(error_stats = summary_vec, 
             colLosses = colLoss_list, 
             fullEsts = fullEsts)
  return(ans)
}

#' @title Validation using interval imputations
#' @param data Matrix or data frame
#' @param rank Rank of decomposition
#' @param leaveout Percent of observations to remove and validate 
#' @param pred.level Confidence of prediction interval
#' @param per.column If not NaN, compute only for that column
xpca_iv = function(data, 
                   rank  = 1,
                    leaveout = 0.25,
                    pred.level=0.90,
                    per.column=NaN){
    # Num cols, rows
    m = nrow(data)
    n = ncol(data)

    # Get obs
    if (!is.na(per.column)) {
        # Get the rows in which that column has no missing
        obs = which(!is.na(data[,per.column])) 
        # Convert to indices
        obs = obs + ((per.column-1)*m)
    }
    else {
        obs = which(!is.na(data)) }
    # Calc num to leave out
    num.leaveout = round(leaveout * length(obs))
    
    # Determine which to leave out
    removed = obs[sample(1:length(obs), num.leaveout)]
    
    # Create matrix with removed entries
    data.obs = data
    data.obs[removed] = NaN

    # Decompose
    fit.xpca = xpca(data.obs, rank=rank)    
    
    # Get each left out entry
    count = 0
    for (e in removed) {
        i = ((e-1) %% m) + 1
        j = floor((e-1) / m) + 1
        ints = getPredictionInterval(i, j,fit.xpca,predLevel=pred.level)
        if (fit.xpca$meanEsts[e] > ints[1] && fit.xpca$meanEsts[e] < ints[2]) {
            count = count + 1
        }
    }

    # Return T if the ratio of the number of imputations
    # that fall within the interval is over the confidence
    ans = list()
    ans$totrust = ((count/num.leaveout) > pred.level)
    ans$within.interval = count
    ans$inds.leftout = removed
    ans$xpca.fit = fit.xpca

    return(ans)
}
