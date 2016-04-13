\name{DESP_OLS_B}
\alias{DESP_OLS_B}
\title{Estimation of B by ordinary least squares}
\description{
  This function estimates the coefficient matrix B column-by-column. Only nonzero coefficients are estimated, using the sparsity pattern provided by the squared partial correlations matrix. gives 
}
\usage{
  DESP_OLS_B(X, SPC)
}
\arguments{
  \item{X}{The design matrix.}
  \item{SPC}{The squared partial correlations matrix.}
}
%\details{
%}
\value{
  The coefficient matrix B.
}

%\references{}

\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link{DESP_SqPartCorr}}}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
