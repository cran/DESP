\name{DESP_SimpleOut}
\alias{DESP_SimpleOut}
\title{Detection of simple outliers}
\description{
  This function aims to detect outliers characterized by a very large Euclidean norm.
}
\usage{
  DESP_SimpleOut(ve, method='MAD', iqr.mult=1.5, mad.mult=2.5, mad.constant=1.4826)
}
\arguments{
  \item{ve}{Vector of Euclidean norms.}
  \item{method}{The method to be used : either based on interquantile range (IQR), or median absolute deviation (MAD, chosen by default).}
  \item{iqr.mult}{The multiplicative factor associated to IQR.}
  \item{mad.mult}{The multiplicative factor associated to MAD.}
  \item{mad.constant}{The scale constant associated to MAD.}
}
\details{
  Denoting the first quartile by Q1, the third one by Q3 and the median by Q2, if the method 'IQR' is chosen, any entry of \code{ve} above Q3 + iqr.mult * (Q3-Q1) is considered as an outlier. If the method 'MAD' is chosen, any entry of \code{ve} above Q2 + mad.mult * MAD is considered as an outlier.
}
\value{
  This function returns the positions of the observations detected as outliers in \code{ve}.
}
%\references{}
\author{Arnak Dalalyan and Samuel Balmand.}

%\seealso{
%}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
