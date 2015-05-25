# DESP/R/DESP_SRL_B.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_SRL_B <-
function(X,lambda,solver="SCS") {
  # estimation of B
  
  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];               # n is the sample size
  p = D[2];               # p is the dimension

  # initialize the matrix B
  B = matrix(,p,p);
  
  # we compute an estimator of each column by the square-root Lasso
  for (j in 1:p)
    {
    beta  = sqR_Lasso(X[,-j],X[,j],lambda,solver);
    B[, j] = append(-beta,1,after=j-1);
    }

  return(B)
}
