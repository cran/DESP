\name{DESP_SRL_B}
\alias{DESP_SRL_B}
\title{Estimation of the coefficient matrix}
\description{
  This function computes the coefficient matrix B estimate, column-by-column.
}
\usage{
  DESP_SRL_B(X, lambda, solver)
}
\arguments{
  \item{X}{The design matrix.}
  \item{lambda}{The penalization parameter.}
  \item{solver}{The solver, SCS by default.}
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
