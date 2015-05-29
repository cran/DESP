\name{sqR_Lasso}
\alias{sqR_Lasso}
\title{
computation of beta that minimize |Y-X*beta|_2 + lambda |beta|_1 (square-root Lasso)
}
\description{
  This function estimates the vector of regression coefficient under sparsity constraints, by square-root Lasso. That is, beta which minimize 
    |Y-X*beta|_2 + lambda |beta|_1.
}
\usage{
  sqR_Lasso(X, Y, lambda, solver)
}
\arguments{
  \item{X}{The matrix of explanatory variables.}
  \item{Y}{The response variable.}
  \item{lambda}{The penalization parameter.}
  \item{solver}{The solver.}
}
\details{
  This method can use the Mosek solver, the Gurobi solver or (by default) the SCS solver.
}
\value{
  The coefficient vector. 
}
%\references{}
\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link[Rmosek]{mosek}},\pkg{gurobi},\code{\link{scsSOCP}}}

\examples{
## set the design matrix
X <- matrix(c(1,0,2,2,1,0,-1,1,1,2,0,1),4,3,byrow=TRUE)
## set the vector of observations
Y <- c(1,0,2,1)
## set the penalty level
lambda <- 1
## compute the Lasso estimate using SCS
sqR_Lasso(X, Y, lambda, solver="SCS")
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
