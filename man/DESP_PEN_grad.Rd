\name{DESP_PEN_grad}
\alias{DESP_PEN_grad}
\title{Steepest descent algorithm for penalized maximum likelihood estimation
}
\description{
  This function implements the steepest descent algorithm with adaptative stepsize and scaled descent direction to solve the maximum likelihood optimization problem and get the inverse of the diagonal of the precision matrix.
}
\usage{
  DESP_PEN_grad(S, B, init, kappa, thresh, stepsize, tol)
}
\arguments{
  \item{S}{The sample covariance matrix.}
  \item{B}{The coefficient matrix.}
  \item{init}{The starting vector of the iteration.}
  \item{kappa}{The tunning paramater.}
  \item{thresh}{The threshold level.}
  \item{stepsize}{The initial step-size.}
  \item{tol}{The gradient magnitude tolerance.}
}
%\details{
%}
\value{
  This function returns the inverse of the diagonal of the precision matrix associated with the sample covariance matrix S as a vector.
}
%\references{}
\author{Arnak Dalalyan and Samuel Balmand.}

%\seealso{
%}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
