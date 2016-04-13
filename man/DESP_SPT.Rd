\name{DESP_SPT}
\alias{DESP_SPT}
\title{Estimation of DESP using shortest path trees}
\description{
  This function estimates the diagonal of the precision matrix using shortest path trees, when the true value of the coefficient matrix B is known or has already been estimated.
}
\usage{
  DESP_SPT(X, B, method = 1)
}
\arguments{
  \item{X}{The design matrix.}
  \item{B}{The coefficient matrix.}
  \item{method}{The method to select the shortest path trees to be considered.}
}
\details{
  The shortest path trees for each connected component could be build under one of the following specifications :
\describe{
  \item{1}{only the presence or absence of edges is considered to build the shortest path trees,}
  \item{2.1}{the root of each shortest path tree is chosen as the node of maximal degree among those in the connected component,}
  \item{2.2}{the root of each shortest path tree - of maximum height equal to 1 - is chosen as the node of maximal degree,}
  \item{3}{for each connected component, the selected shortest path tree is the maximum weighted tree among all shortest path trees.}
}
}
\value{
  This function returns the diagonal of the precision matrix associated with X as a vector.
}
%\references{}
\author{Arnak Dalalyan and Samuel Balmand.}

%\seealso{
%}

\examples{
## build the true precision matrix
p <- 12 # number of variables
Omega <- 2*diag(p)
Omega[1,1] <- p 
Omega[1,2:p] <- 2/sqrt(2)
Omega[2:p,1] <- 2/sqrt(2)
## compute the true diagonal of the precision matrix
Phi <- 1/diag(Omega)
## generate the design matrix from a zero-mean Gaussian distribution
require(MASS)
n <- 100 # sample size
X <- mvrnorm(n,rep.int(0,p),ginv(Omega))
## compute the sample mean
barX <-  apply(X,2,mean)
## subtract the mean from all the rows
X <- t(t(X)-barX)
## estimate the coefficient matrix 
B <- DESP_SRL_B(X,lambda=sqrt(2*log(p)))
## compute the squared partial correlations
SPC <- DESP_SqPartCorr(B,n)
## re-estimate the coefficient matrix by ordinary least squares
B_OLS <- DESP_OLS_B(X,SPC)
## estimate the diagonal of the precision matrix and get its inverse
hatPhiSPT <- 1/DESP_SPT(X,B_OLS,method=2.1)
## measure the performance of the estimation using l2 vector norm
sqrt(sum((Phi-hatPhiSPT)^2))
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
