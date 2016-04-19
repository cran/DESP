\name{desp.cv}
\alias{desp.cv}
\title{Selection of the tuning parameters of desp by v-fold cross-validation}
\description{
  This function returns the precision matrix and the expectation associated to the data matrix X using \code{\link{desp}} and choosing the tuning parameters \eqn{\lambda}{lambda} and \eqn{\gamma}{gamma} by \code{v}-fold cross-validation that uses a robust loss function. The expression of the loss function is provided in the companion vignette.
}
\usage{
  desp.cv(X, v, lambda.range, gamma.range, settings=NULL)
}
\arguments{
  \item{X}{The data matrix.}
  \item{v}{The number of folds.}
  \item{lambda.range}{The range of the penalization parameter \eqn{\lambda}{lambda} that encourages robustness.}
  \item{gamma.range}{The range of the penalization parameter \eqn{\gamma}{gamma} that promotes sparsity.}
  \item{settings}{A list including all the parameters needed for the estimation. Please refer to the documentation of the function \code{\link{desp}} to get more details.}
}
%\details{
%}
\value{
  \code{desp.cv} returns an object with S3 class "desp.cv" containing the estimated parameters along with the selected values of the tuning parameters, with components:
  \item{Omega}{The precision matrix.}
  \item{mu}{The expectation vector.}
  \item{Theta}{The matrix corresponding to outliers.}
  \item{lambda}{The selected lambda.}
  \item{gamma}{The selected gamma.}
}
%\references{
%}
\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link{desp}}}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
