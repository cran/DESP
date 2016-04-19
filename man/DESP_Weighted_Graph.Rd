\name{DESP_Weighted_Graph}
\alias{DESP_Weighted_Graph}
\title{Graph representation from the matrix B}
\description{
  This function build a weighted undirected graph without loops from the matrix \code{B}. The weights of the edges depend on squared partial correlations.
}
\usage{
  DESP_Weighted_Graph(B, n)
}
\arguments{
  \item{B}{The coefficient matrix.}
  \item{n}{The sample size.}
}
%\details{
%}
\value{
  A new undirected graphNEL object.
}

%\references{}

\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link[graph]{sparseM2Graph}}}

%\examples{
%}


%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
