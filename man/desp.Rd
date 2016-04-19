\name{desp}
\alias{desp}
\title{Robust estimation of a sparse precision matrix}
\description{
  This function computes robust estimates of the precision matrix \eqn{\Omega}{Omega} and of the expectation \eqn{\mu}{mu} corresponding to the data matrix \code{X} by solving a convex optimization problem. It implements the algorithms proposed by Balmand and Dalalyan (2015a, 2015b).
}
\usage{
  desp(X,lambda=0,gamma=0,settings=NULL)
}
\arguments{
  \item{X}{The data matrix.}
  \item{lambda}{The penalization parameter that encourages robustness.}
  \item{gamma}{The penalization parameter that promotes sparsity.}
  \item{settings}{A list including all the parameters needed for the estimation.
    \describe{
      \item{\code{nThreads}}{number of threads for parallel processing, default 1}
      \item{\code{thres}}{threshold used as a stopping criterion, the computation stops when the variation of the objective function is below \code{thres}, default 1e-6}
      \item{\code{maxIter}}{maximum number of iterations, default 2000}
      \item{\code{diagElem}}{estimator of the diagonal entries, default 'RV' (residual variance). The other available estimators}
      \item{\code{OLS}}{whether to re-estimate the coefficients of regression by ordinary least squares on selected variables only or not, default TRUE}
      \item{\code{OLS.min}}{the number of observations detected as inliers below which the coefficients of regression are not re-estimated, used only when \code{OLS} is TRUE, half of the sample size by default}
      \item{\code{refit}}{whether to re-estimate the coefficients of regression after having removed the detected outliers, default TRUE}
      \item{\code{refit.min}}{the number of observations detected as inliers below which the coefficients of regression are not re-estimated, used only when \code{refit} is TRUE, half of the sample size by default}
      \item{\code{sqrLasso}}{method used to solve the square-root Lasso problem, default 'CD'}
      \item{\code{sto}}{whether a randomized algorithm (stochastic coordinate descent) have to be used when choosing the coordinate descent method. By default, this parameter is set to '0', that means that the coordinates are updated in the order in which the corresponding variables appear in X. Another option would be '2', the coordinates are all updated but in a uniformly random order. The last option (experimental) would be '1', in this case the sole coordinate to be updated is chosen uniformly at random at each iteration.}
      \item{\code{outParCorrB}}{whether the partial correlation matrix and the matrix of the coefficients of regression must be outputted or not, default FALSE}
      \item{\code{posDef}}{whether to ensure that the precision matrix is positive definite or not, default TRUE}
      \item{\code{SML.method}}{when \code{diagElem} is 'SML', the method to choose the path between two vertices, default 'MST'}
      \item{\code{SPT.spec}}{when \code{diagElem} is 'SML' and when \code{SML.method} is 'SPT', the specifications to build the shortest path trees for each connected component, default 1}
      \item{}{parameters to be used when \code{diagElem} is 'PML' :
        \describe{
        \item{\code{PML.thresh}}{threshold level, default 0.05}
        \item{\code{PML.kappa}}{tunning paramater, default 0.05}
        \item{\code{PML.tol}}{gradient magnitude tolerance, default 1e-5}
        }
      }
      \item{\code{rmSO}}{whether the outliers characterized by a very large Euclidean norm have to be excluded first, default TRUE}
      \item{\code{rmSO.method}}{method used to detect rough outliers: either based on interquantile range (IQR), or median absolute deviation, used only when \code{rmSO} is TRUE, default 'MAD'}
      \item{\code{rmSO.mad.mult}}{multiplicative factor associated to MAD, default 5}
      \item{\code{rmSO.mad.constant}}{scale constant associated to MAD, default 1.4826}
      \item{\code{rmSO.iqr.mult}}{multiplicative factor associated to IQR, default 3}
    }
  }
}
%\details{
%}
\value{
  \code{desp} returns an object with S3 class "desp" containing the estimated parameters, with components:
  \item{Omega}{The precision matrix.}
  \item{mu}{The expectation vector.}
  \item{Theta}{The matrix corresponding to outliers.}
  \item{Psi}{The matrix of partial correlations, optional.}
  \item{B}{The matrix of the coefficients of regression, optional.}
  \item{simpleOut}{The outliers characterized by a very large Euclidean norm, optional. This observations have been excluded and do not appear in \code{Theta}}
}
\references{
  Balmand, S. and Dalalyan, A. S. (2015a): On estimation of the diagonal elements of a sparse precision matrix, Posted as http://arxiv.org/abs/1504.04696.

  Balmand, S. and Dalalyan, A. S. (2015b): Convex programming approach to robust estimation of a multivariate Gaussian model, Posted as http://arxiv.org/abs/1512.04734.
}
\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link{sqR_Lasso}}}

\examples{
## set the dimension and the sample size
p <- 50
n <- 40
## define the true precision matrix
Omega <- 2*diag(p)
Omega[1,1] <- p 
Omega[1,2:p] <- 2/sqrt(2)
Omega[2:p,1] <- 2/sqrt(2)
## generate the distribution
set.seed(1)
X <- MASS::mvrnorm(n, rep.int(0,p), MASS::ginv(Omega))
## define the settings
settings <- list("diagElem"="AD", "OLS"=FALSE)
## estimate the parameters of the distribution
params <- desp(X,lambda=0,gamma=sqrt(2*log(p)),settings=settings)
## error of estimation measured by Frobenius norm
sqrt(sum((Omega - params$Omega)^2))
## increase the sample size and generate the distribution again
n <- 1000
X <- MASS::mvrnorm(n, rep.int(0,p), MASS::ginv(Omega))
## estimate the parameters of the distribution
params <- desp(X,lambda=0,gamma=sqrt(2*log(p)),settings=settings)
## error of estimation measured by Frobenius norm
sqrt(sum((Omega - params$Omega)^2))
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
