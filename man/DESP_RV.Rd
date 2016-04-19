\name{DESP_RV}
\alias{DESP_RV}
\title{Estimation of DESP by residual variance}
\description{
  This function estimates the diagonal of the precision matrix by residual variance when the true value of the coefficient matrix \code{B} is known or has already been estimated. The observations of the data matrix \code{X} are assumed to have zero mean.
}
\usage{
  DESP_RV(X, B, Theta = NULL)
}
\arguments{
  \item{X}{The data matrix.}
  \item{B}{The coefficient matrix.}
  \item{Theta}{The matrix orresponding to outliers.}
}
\details{
  When Theta is not NULL, we consider an additive contamination model. We assume that X = Y + E is observed, denoting the outlier-free data by Y and the matrix of errors by E. In this case, the matrix Theta should be equal to E * B.
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
hatPhiRV <- 1/DESP_RV(X,B_OLS)
## measure the performance of the estimation using l2 vector norm
sqrt(sum((Phi-hatPhiRV)^2))
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
