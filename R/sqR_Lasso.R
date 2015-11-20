# DESP/R/sqR_Lasso.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

sqR_Lasso <-
function(X,Y,lambda,solver="SCS") {
  # computation of beta that minimize |Y-X*beta|_2 + lambda |beta|_1 (square-root Lasso)
  
  if (solver=="SCS"){
    r <- try(.Call("scs_sqR_Lasso",as.matrix(X),as.double(Y),as.double(lambda)), silent = TRUE)
    if (r$status != 1) {
      stop (paste("SCS status :",r$status))
    }
    if ( inherits (r , "try-error")) {
      stop ("SCS failed somehow !")
    }
    beta <- as.matrix(r$beta)
  }
  else{

  # read the sample size and the number of variables
  D <- dim(X)
  n <- D[1]               # n is the sample size
  p <- D[2]               # p is the dimension

  # normalize the columns of X by dividing them by their norm
  vn <- 1/sqrt(apply(X^2,2,mean))
  X_ <- tcrossprod(X,diag(vn))
  
  # we now form the optimization problem, which takes the following form. 
  if (solver=="Mosek"){
    # The variables are xx=c(z=abs(beta),beta,v=Y-X*beta,t=norm(v))
    # 
  
    # the cost function
    qo <- list(sense = "min")
    qo$c <- c((1:p)*0+lambda, (1:(p+n))*0, 1)
  
    # these lines introduce the constraints |beta_{j}| = z_{j}
    A <- rBind(cBind( Diagonal(p), Diagonal(p), Matrix(0, p, n+1)), cBind( Diagonal(p) , -Diagonal(p), Matrix(0, p, n+1)))
    blc <- (1:(2*p))*0
    buc <- (1:(2*p))*0+Inf
  
    # these are the constraints v+X*beta=Y
  
    A <- rBind(A, cBind(Matrix(0,n,p), X_, Diagonal(n), Matrix(0,n,1)))
    blc <- rbind(t(t(blc)), t(t(Y)));
    buc <- rbind(t(t(buc)), t(t(Y)));
  
    # we finally define the quadratic constraint |v|_2 <= t.
  
    qo$cones <- list ("QUAD", c(p+p+n+1 , 2*p+1:n))
  
    qo$A <- A;
    qo$bc <- rbind(t(blc),t(buc));
    qo$bx <- rbind( (1:(2*p+n+1))*0-Inf, (1:(2*p+n+1))*0+Inf);
  
    NUMCONES <- 1
    qo$cones <- matrix( list(), nrow =2, ncol = NUMCONES )
    rownames( qo$cones ) = c("type","sub")
    qo$cones[ ,1] = list("QUAD", c(2*p+n+1 , 2*p+1:n))
  
    r <- try(Rmosek::mosek(qo,list( verbose =0 )), silent = TRUE )
    if (r$response$code == 1001) {
      stop (paste("Rmosek :",r$response$msg))
    }
    if ( inherits (r , "try-error")) {
      stop ("Rmosek failed somehow !")
    }
    x <- r$sol$itr$xx
  }
  if (solver=="Gurobi"){

    # A matrix
    # these line introduce the constraints |beta_{j}| = z_{j}
    A <- rBind(cBind( Diagonal(p), Diagonal(p), Matrix(0, p, n+1)), cBind( Diagonal(p) , -Diagonal(p), Matrix(0, p, n+1)))
    # these are the constraints v+X*beta=Y
    A <- rBind(A, cBind(Matrix(0,n,p), X_, Diagonal(n), Matrix(0,n,1)))
    
    params <- list(OutputFlag=0)#LogToConsole=0, LogFile=""
    
    model <- list()
    model$A      <- A
    model$cones  <- list(as.list(c(p+p+n+1,2*p+1:n)))
    model$obj    <- c((1:p)*0+lambda, (1:(p+n))*0, 1) # c
    model$rhs    <- c((1:p)*0,(1:p)*0,Y)
    model$lb     <- (1:(2*p+n+1))*0-Inf
    #model$ub     <- (1:(2*p+n+1))*0+Inf # default 
    model$sense  <- c(rep('>=',p), rep('>=',p), rep('=',n))
    #model$modelsense <- "min" # default
    #model$vtype      <- 'C' # default
    
    r <- try(gurobi::gurobi(model, params), silent = TRUE)
    if (r$status != "OPTIMAL") {
      stop (paste("Gurobi :",r$status))
    }
    if ( inherits (r , "try-error")) {
      stop ("Gurobi failed somehow !")
    }

    x <- r$x
  }
  beta <- as.matrix(x[p+1:p] * vn)
  }

  return(beta)

}
