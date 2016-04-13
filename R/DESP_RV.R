# DESP/R/DESP_RV.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_RV <-
function(X,B) {
  # estimation of the diagonal of the precision matrix by residual variance when the true value of B is known or has already been estimated
  
  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];               # n is the sample size
  p = D[2];               # p is the dimension

  # initialize Phi
  Phi = rep(0,p);

  for (j in 1:p)
	{
	beta  = -B[-j,j];
	sigma2 = mean((X[,j]-X[,-j]%*%beta)^2);
	Phi[j]  = ifelse(sigma2<1,sigma2,1); 
	}

  return(1/Phi);
}
