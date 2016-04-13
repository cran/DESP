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

\examples{
## This example is the same that the one given for sqR_Lasso() function, 
## with the difference that the optimization problem in R is set up in R,
## before calling scsSOCP() function.
##
## set the design matrix
X <- matrix(c(1,0,2,2,1,0,-1,1,1,2,0,1),4,3,byrow=TRUE)
## set the vector of observations
Y <- c(1,0,2,1)
## set the penalty level
lambda <- 1
## compute the square-root Lasso estimate using SCS
## computation of beta that minimize |Y-X*beta|_2 + lambda |beta|_1
##
## read the sample size and the number of variables
D <- dim(X)
n <- D[1]               # n is the sample size
p <- D[2]               # p is the dimension
## normalize the columns of X by dividing them by their norm
vn <- 1/sqrt(apply(X^2,2,mean))
X_ <- tcrossprod(X,diag(vn))
## minimize |Y-X*beta|_2 + lambda |beta|_1 (square-root Lasso)
## using the variables : x=c(z=abs(beta),beta,v=Y-X*beta,t=norm(v))
## the optimization problem could be rewritten :
##  minimize
##        c^T x
##  subject to
##        v + X beta  = Y    (1)
##        z + beta  >= 0     (2)
##        z - beta  >= 0     (3)
##        |v|_2^ 2 <= t ^ 2  (4) (a second-order cone constraint)
##
## to solve this problem using SCS, we write
##  minimize
##        c^T x
##  subject to
##        b - A x  = s
##        s in K
## where K is a product of zero, linear and second-order cones, 
## containing the following informations (in this order) :
## f = n : number of linear equality constraints (cf. (1))
## l = 2*p : length of LP cone (cf. (2) and (3))
## *q = [1+n] : array of second-order cone constraints (cf. (4))
## qsize = 1 : length of SOC array (number of second-order cones)
##
require(Matrix)
## A matrix
## these lines introduce the constraints (1) : v+X*beta=Y
A <- cBind(Matrix(0,n,p), X_, Diagonal(n), Matrix(0,n,1))
## these lines introduce the constraints (2) and (3) : |beta_{j}| = z_{j}
A <- rBind(A,cBind( -Diagonal(p), -Diagonal(p), Matrix(0, p, n+1)), 
cBind( -Diagonal(p) , Diagonal(p), Matrix(0, p, n+1)))
## these lines introduce the quadratic constraint (4) : |v|_2 <= t
A <- rBind(A, cBind(Matrix(0, 1, 2*p+n),-1), cBind(Matrix(0, n, 2*p), Diagonal(n), Matrix(0, n, 1)))
##
## b vector
b <- c(Y,(1:(2*p+1+n))*0)
##
## c vector : the cost function
c <- c((1:p)*0+lambda, (1:(p+n))*0, 1)
##
## K vector of second-order cone constraints
K <- c(1+n)
##
model <- list()
model$Ax <- A@x # access to the slot containing the non-zero elements of the matrix
model$Ai <- A@i # access to the slot containing the row index of class dgCMatrix,
		# subclass of CsparseMatrix
model$Ap <- A@p # access to the slot containing the column pointer of class dgCMatrix,
		# subclass of CsparseMatrix
model$Am <- 2*n+2*p+1 # A number of rows
model$An <- 2*p+n+1 # A number of columns
model$b <- b
model$c <- c
model$Kf <- n
model$Kl <- 2*p
model$Kq <- K
model$Kqsize <- 1
model$n <- n
model$p <- p
model$settings <- list(normalize=0,verbose=0)
##
r <- try(scsSOCP(model), silent = TRUE)
if (r$status != 1) {
  stop (paste("SCS status :",r$status))
}
if ( inherits (r , "try-error")) {
  stop ("SCS failed somehow !")
}
## get beta, the vector of the coefficients of regression
as.matrix(r$x[p+1:p] * vn)
}

\keyword{ SOCP }
%\keyword{ ~kwd2 }




