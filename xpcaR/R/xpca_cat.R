#' @import data.table mltools
NULL

softMax_fix <- function(xpca_fit, orgColNames){
  ests = xpca_fit$meanEsts
  expandedColNames = colnames(ests)
  if(length(expandedColNames) == length(orgColNames)){
    warning("length of orgColNames matches number of columns in fit. No softMax_fit performed")
    return(xpca_fit)
  }
  for(j in seq_along(orgColNames)){
    this_name = orgColNames[j]
    is_this_var = which(grepl(this_name, expandedColNames) )
    if(length(is_this_var) > 1){
      rs <- rowSums(ests[,is_this_var])
      denom = NULL
      for(j in 1:length(is_this_var))
        denom = cbind(denom, rs)
      ests[,is_this_var] <- ests[,is_this_var] / denom
    }
  }
  xpca_fit$meanEsts_softMax <- ests
  return(xpca_fit)
}

expandData <- function(df){
  dt <- data.table::as.data.table(df)
  ans = mltools::one_hot(dt, dropUnusedLevels = T)
  ans = as.data.frame(ans)
  return(ans)
}

xpca_cv_cat <- function(data, 
                   folds = 4, 
                   methods = c("pca", "xpca", "coca"), 
                   rank  = 1, 
                   lossFun = cv_rmse){
  colnames(data) <- paste0("_", colnames(data))
  org_names = colnames(data)
  data = expandData(data)
  cnames = colnames(data)
  # Determining which data is observed
  data_mat = as.matrix(data)
  is_obs = which(!is.na(data_mat))
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
  xpca_cat_ests = data_mat + NA
  coca_ests = data_mat + NA
  coca_cat_ests = data_mat + NA
  pca_ests  = data_mat + NA
  
  # Partitioning data into folds, computing decompositions 
  # and filling in ests matrices
  for(i in 1:folds){
    this_data = data_mat
    these_indices = partitions[[i]]
    this_data[these_indices] = NA
    if(any(methods == "xpca")){
      fit = xpca(this_data, rank = rank)
      fit = softMax_fix(fit, org_names)
      xpca_ests[these_indices] = fit$meanEsts[these_indices]
      xpca_cat_ests[these_indices] = fit$meanEsts_softMax[these_indices]
    }
    if(any(methods == "pca")){
      fit = pca(this_data, rank = rank)
      pca_ests[these_indices] = fit$meanEsts[these_indices]
    }
    if(any(methods == "coca")){
      fit = coca(this_data, rank = rank)
      fit = softMax_fix(fit, org_names)
      coca_ests[these_indices] = fit$meanEsts[these_indices]
      coca_cat_ests[these_indices] = fit$meanEsts_softMax[these_indices]
    }
  }
  # Creating summary stat for each method and saving full ests
  summary_vec = NULL
  summary_vec_names = NULL
  full_ests = list()
  colLoss_list = list()
  
  if(any(methods == "xpca")){
    loss_vec = NULL
    loss_vec_cat = NULL
    for(j in 1:ncol(data_mat)){ 
      loss_vec[j] = lossFun( data_mat[,j], xpca_ests[,j], j)
      loss_vec_cat[j] = lossFun( data_mat[,j], xpca_cat_ests[,j], j)
    }
    names(loss_vec) <- cnames
    names(loss_vec_cat) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec), mean(loss_vec_cat))
    summary_vec_names = c(summary_vec_names, "xpca", "xpca_cat")
    colLoss_list[["xpca"]] = loss_vec
    colLoss_list[["xpca_cat"]] = loss_vec_cat
    full_ests[["xpca"]] = xpca_ests
    full_ests[["xpca_cat"]] = xpca_cat_ests
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
    loss_vec_cat = NULL
    for(j in 1:ncol(data_mat)){ 
      loss_vec[j] = lossFun( data_mat[,j], coca_ests[,j], j)
      loss_vec_cat[j] = lossFun( data_mat[,j], coca_cat_ests[,j], j)
    }
    names(loss_vec) <- cnames
    names(loss_vec_cat) <- cnames
    summary_vec = c(summary_vec, mean(loss_vec), mean(loss_vec_cat))
    summary_vec_names = c(summary_vec_names, "coca", "coca_cat")
    colLoss_list[["coca"]] = loss_vec
    colLoss_list[["coca_cat"]] = loss_vec_cat
    full_ests[["coca"]] = xpca_ests
    full_ests[["coca_cat"]] = xpca_cat_ests
  }
  # Prepping output
  colLossMat = NULL
  for(i in seq_along(colLoss_list)){
    colLossMat <- cbind(colLossMat, colLoss_list[[i]])
  }
  colnames(colLossMat) = names(colLoss_list)
  names(summary_vec) = summary_vec_names
  ans = list(error_stats = summary_vec, 
             colLosses = colLossMat,
             fullEsts = full_ests)
  return(ans)
}

dropRareOutcome_one <- function(factor, min = 10){
  tabledValues <- table(factor)
  isTooSmall <- which(tabledValues < min)
  for(i in seq_along(isTooSmall)){
    this_ind = isTooSmall[i]
    this_val = names(tabledValues)[this_ind]
    factor[factor == this_val] = NA
  }
  return(factor)
}

dropRareOutcome <- function(data, min = 10){
  for(j in 1:ncol(data)){
    if(is.factor(data[,j])){
      data[,j] = dropRareOutcome_one(data[,j], min = min)
    }
  }
  return(data)
}
