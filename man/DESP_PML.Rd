\name{DESP_PML}
\alias{DESP_PML}
\title{Estimation of DESP by penalized likelihood minimization}
\description{
  This function estimates the diagonal of the precision matrix by penalized likelihood minimization, when the true value of the coefficient matrix B is known or has already been estimated.
}
\usage{
  DESP_PML(X, B, thresh, kappa, tol)
}
\arguments{
  \item{X}{The design matrix.}
  \item{B}{The coefficient matrix.}
  \item{thresh}{The threshold level.}
  \item{kappa}{The tunning paramater.}
  \item{tol}{The gradient magnitude tolerance.}
}
%\details{
%}
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
hatPhiPML <- 1/DESP_PML(X,B,thresh=0.05,kappa=0.1, tol=1e-5)
## measure the performance of the estimation using l2 vector norm
sqrt(sum((Phi-hatPhiPML)^2))
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
