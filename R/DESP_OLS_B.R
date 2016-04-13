# DESP/R/DESP_OLS.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_OLS_B <-
function(X,SPC) {
  # estimation of B by ordinary least squares
  
  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];           # n is the sample size
  p = D[2];           # p is the dimension

  # compute the sample covariance matrix
  S = crossprod(X)/n;

  # initialize the matrix B
  B = matrix(,p,p);
  
  # we compute an estimator of each line for every selected (non-zero) coefficient by OLS
  for (j in 1:p)
    {
    SPCc  = as.matrix(SPC)[-j,j];
    select = which(SPCc != 0);
    beta = c(1:(p-1))*0;
    if(length(select)!=0){
      beta[select] = crossprod(ginv(crossprod((X[,-j])[,select])) , crossprod((X[,-j])[,select] , X[,j]));
    }

    B[, j] = append(-beta,1,after=j-1);
    }

  return(B)
}
