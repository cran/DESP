\name{scsSOCP}
\alias{scsSOCP}
\title{
solve a second-order cone program using SCS
}
\description{
  This function solves a second-order cone program using the C package SCS that solves convex cone problems via operator splitting. We use the direct method of the C package SCS that also proposes an indirect method.
  This method solves an optimization problem of the form :
    minimize
          c^T x
    subject to
          b - A x  = s,
          s in K,
    where K is a product of zero, linear and second-order cones, containing the following informations (in this order) :
    f : number of linear equality constraints,
    l : length of LP cone,
    q : array of second-order cone constraints,
    qsize : length of SOC array (number of second-order cones).
}
\usage{
  scsSOCP(model)
}
\arguments{
  \item{model}{The model representing the convex problem to solve is a list containing the components that follow. The sparse matrix A is stored using the compressed sparse column (CSC) format.
    \describe{
    \item{Ax}{non-zero elements of the matrix A, size : nnz A}
    \item{Ai}{row index of non-zero elements of the matrix A, size : nnz A}
    \item{Ap}{column pointer for non-zero elements of the matrix A, size: An + 1}
    \item{Am}{number of rows of the matrix A}
    \item{An}{number of columns of the matrix A}
    \item{b}{dense array for b (size Am)}
    \item{c}{dense array for c (size An)}
    \item{Kf}{number of linear equality constraints}
    \item{Kl}{length of LP cone}
    \item{Kq}{array of second-order cone constraints}
    \item{Kqsize}{length of SOC array (number of second-order cones)}
    \item{settings}{list of settings composed of :
      \describe{
      \item{alpha}{relaxation parameter, default 1.8}
      \item{rho_x}{x equality constraint scaling, default 1e-3}
      \item{max_iters}{maximum iterations to take, default 2500}
      \item{eps}{convergence tolerance, default 1e-3}
      \item{cg_rate}{for indirect, tolerance goes down like (1/iter)^cg_rate, default 2}
      \item{verbose}{boolean, write out progress, default 1(true)}
      \item{normalize}{boolean, heuristic data rescaling, default 1(true)}
      \item{scale}{if normalized, rescales by this factor, default 5}
      \item{warm_start}{boolean, warm start (put initial guess in Sol struct), default 0(false)}
      }
    }
    }
  }
}
\details{
  More documentation on SCS can be found at \url{https://github.com/cvxgrp/scs}.
}
\value{
  \code{scsSOCP} returns a list containing the optimal solution, with components:
  \item{x}{The value of the primal variable for the best solution.}
  \item{status}{The status of the optimization, returned as an integer : 1 in case of a solved program. A negative status corresponds to an infeasible, unbounded, inaccurate, interrupted, failed or indeterminate problem. Refer to the SCS development pages for more documentation on status codes.}
}

\references{
  O'Donoghue, B. and Chu, E. and Parikh, N. and Boyd, S. (2015): Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding.
}

\author{Arnak Dalalyan and Samuel Balmand.}

%\seealso{
%}

%\examples{
%}

\keyword{ SOCP }
%\keyword{ ~kwd2 }




