\name{DESP_SRL_B}
\alias{DESP_SRL_B}
\title{Estimation of the coefficient matrix}
\description{
  This function computes the coefficient matrix \code{B} estimate, column-by-column. The observations of the data matrix \code{X} are assumed to have zero mean.
}
\usage{
  DESP_SRL_B(X, lambda, solver = 'CD', sto = FALSE, nThreads = 1)
}
\arguments{
  \item{X}{The data matrix.}
  \item{lambda}{The penalization parameter that promotes sparsity.}
  \item{solver}{The solver.
    A string indicating the solver to use.
    \describe{
      \item{\code{"CD"}}{specifies the coordinate descent algorithm;}
      \item{\code{"SCS"}}{specifies the free Splitting Conic Solver;}
      \item{\code{"Gurobi"}}{specifies the commercial Gurobi solver;}
      \item{\code{"Mosek"}}{specifies the commercial Mosek solver.}
    }
    The default is \code{"CD"}.}
  \item{sto}{Indicates whether a randomized algorithm (stochastic coordinate descent) have to be used when choosing the coordinate descent method, default FALSE.}
  \item{nThreads}{Set the number of threads for parallel processing when choosing the coordinate descent method or the SCS solver, default 1.}
}
%\details{
%}
\value{
  The coefficient matrix.
}
%\references{}
\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link{sqR_Lasso}}}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
