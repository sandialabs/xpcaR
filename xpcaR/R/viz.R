vis1factor <- function(fact_ind, factors, 
                       use_rnames = T,
                       ...){
  these_vals = factors[,fact_ind]
  argList = list(...)
  argList$col = 0
  argList$x = these_vals
  do.call(plot, argList)
  rnames = rownames(factors)
  if(use_rnames & !is.null(rnames))
    text(1:length(these_vals), these_vals, 
         rnames, 
         cex = 0.5, ...)
  else
    points(these_vals, pch = 16, ...)
  lines(c(0, length(these_vals) + 1), c(0,0), 
        col = 'red', lwd = 2)
}

#' @title Visualizing Decomposition Factors
#' @param decomp An xpcaR::pca, coca or XPCA object
#' @param usenames Should row/column names be plotted?
#' @param ... additional arguments passed to `plot`
#' @export
visFactors = function(decomp, usenames = T, ...){
  row_fctrs = decomp$A
  col_fctrs = decomp$B

  nFctrs = ncol(col_fctrs)
  par(mfrow = c(nFctrs, 2))
  for(i in seq_len(nFctrs)){
    fctrName = paste("PC", i)
    vis1factor(i, col_fctrs, 
               ylab = fctrName, 
               main = paste0(fctrName, ": Columns"),
               use_rnames = usenames, 
               ...)
    vis1factor(i, row_fctrs, 
               ylab = fctrName, 
               main = paste0(fctrName, ": Rows"),
               use_rnames = usenames,
               ...)
  }
}

