% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xPCA.R
\name{xpca}
\alias{xpca}
\title{XPCA}
\usage{
xpca(
  data,
  rank,
  tol = 0.1,
  maxIters = 1000,
  opt = "AN",
  reg_val = 0.5,
  post_svd = T,
  gridSize = nrow(data)/10
)
}
\arguments{
\item{data}{Data matrix or data.frame}

\item{rank}{Rank of decomposition}

\item{tol}{Tolerance; if change in LLK is below this value, algorithm terminates}

\item{maxIters}{Maximum iterations}

\item{opt}{Optimization algorithm to use. Options are "LBFGS" or "AN" (Alternating Newton's)}

\item{reg_val}{Value used to regularize max intervals. Prevents Hauck-Donner effect.}

\item{post_svd}{Should svd be used to standardize decomposition?}

\item{gridSize}{Number of points used to approximate mean function}
}
\value{
A list with the following fields:

   - \code{A}: Low-rank representation of rows
   
   - \code{B}: Low-rank representation of columns
   
   - \code{fittedEsts}: Data estimates of low-rank 
   representation of full data matrix
   
Other fields are for internal use.
}
\description{
Computes low-rank XPCA decomposition of data. 
Can be used with missing data. If all data is continuous, 
estimates will be very close to those provided by COCA. 
If marginal distribution of all columns is approximately 
Gaussian, both COCA + XPCA estimates will be very close 
to those provided by PCA.
}
\details{
Copula-based decomposition, 
where all data is treated as ordinal. 

For technical details, see:

C. Anderson-Bergman, T. G. Kolda, K. Kincher-Winoto. 
XPCA: Extending PCA for a 
Combination of Discrete and Continuous Variables. 
arXiv:1808.07510, 2018
}
\examples{
data = simProblem()$data
fit = xpca(data, rank = 2)
}
