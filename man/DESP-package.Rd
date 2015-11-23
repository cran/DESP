\name{DESP-package}
\alias{DESP-package}
\alias{DESP}
\docType{package}
\title{Estimation of Diagonal Elements of Sparse Precision-Matrices}
\description{
  This package provides several methods of estimation of the diagonal elements of a sparse precision matrix. It also includes functions for estimating the off-diagonal elements of the latter by the square-root Lasso. 
}
\details{
\tabular{ll}{
Package: \tab DESP\cr
Type: \tab Package\cr
Version: \tab 0.1-6\cr
Date: \tab 2015-11-23\cr
License: \tab GPL-3\cr
}
\describe{
\item{\code{DESP_RV}}{Estimation of DESP by residual variance}
\item{\code{DESP_RML}}{Estimation of DESP by likelihood maximization}
\item{\code{DESP_SPT}}{Estimation of DESP using shortest path trees}
\item{\code{DESP_MST}}{Estimation of DESP using minimum spanning trees}
\item{\code{DESP_PML}}{Estimation of DESP by penalized likelihood minimization}
\item{\code{sqR_Lasso}}{Solve the square-root LASSO}
\item{\code{scsSOCP}}{Solve a second-order cone program using SCS}
}
}
\author{
Arnak Dalalyan and Samuel Balmand

Maintainer: Samuel Balmand <Samuel.Balmand@ensg.eu>
}
\references{
  Balmand, S. and Dalalyan, A. S. (2015): On estimation of the diagonal elements of a sparse precision matrix, Posted as http://arxiv.org/abs/1504.04696.
}

\keyword{ package }
\seealso{
\code{\link[Rmosek]{mosek}},\pkg{gurobi},\code{\link{scsSOCP}}, \code{\link[graph]{connComp}}, \code{\link[RBGL]{mstree.kruskal}}, \code{\link[RBGL]{dijkstra.sp}} 
}

