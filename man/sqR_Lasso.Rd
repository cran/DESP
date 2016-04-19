\name{sqR_Lasso}
\alias{sqR_Lasso}
\title{
computation of beta that minimize |Y-X*beta|_2 + lambda |beta|_1 (square-root Lasso)
}
\description{
  This function estimates the vector of regression coefficient under sparsity constraints, by square-root Lasso. That is, \eqn{\beta}{beta} which minimize 
    \deqn{{\|Y-X \beta\|}_2 + \lambda {\|\beta\|}_1.}{%
    |Y-X*beta|_2 + lambda |beta|_1.}
}
\usage{
  sqR_Lasso(X, Y, lambda, solver = 'CD', sto = '0')
}
\arguments{
  \item{X}{The matrix of explanatory variables (must be a double-precision matrix).}
  \item{Y}{The response variable.}
  \item{lambda}{The penalization parameter.}
  \item{solver}{The solver.
    A string indicating the solver to use.
    \describe{
      \item{\code{"CD"}}{specifies the coordinate descent algorithm;}
      \item{\code{"SCS"}}{specifies the free Splitting Conic Solver;}
      \item{\code{"Gurobi"}}{specifies the commercial Gurobi solver;}
      \item{\code{"Mosek"}}{specifies the commercial Mosek solver.}
    }
    The default is \code{"CD"}.}
  \item{sto}{Indicates whether a randomized algorithm (stochastic coordinate descent) have to be used when choosing the coordinate descent method. By default, this parameter is set to '0', that means that the coordinates are updated in the order in which the corresponding variables appear in X. Another option would be '2', the coordinates are all updated but in a uniformly random order. The last option (experimental) would be '1', in this case the sole coordinate to be updated is chosen uniformly at random at each iteration.}
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
## compute the square-root Lasso estimate using SCS
## get beta, the vector of the coefficients of regression
sqR_Lasso(X, Y, lambda, solver="SCS")
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
