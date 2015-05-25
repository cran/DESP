\name{DESP_SqPartCorr}
\alias{DESP_SqPartCorr}
\title{Squared partial correlations computation}
\description{
  This function calculates the squared partial correlations. We keep only the elements upper than a threshold. For example, if thresh=0.01, therefore SPC[i,j] equal 0 means that the proportion of the variance of the variable i not explained by the other variables, which is explained by the variable j is less than 1%.
}
\usage{
  DESP_SqPartCorr(B, n)
}
\arguments{
  \item{B}{The coefficient matrix.}
  \item{n}{The sample size.}
}
%\details{
%}
\value{
  The matrix of squared partial correlations.
}
%\references{}
\author{Arnak Dalalyan and Samuel Balmand.}

%\seealso{
%}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
