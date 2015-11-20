/*
 *  DESP/src/scsSqR_Lasso_solve.h by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (version 3) as published by
 *  the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 */

#ifndef SCSSQR_LASSO
#define SCSSQR_LASSO

#include "scs.h"
#include "linsys/amatrix.h"

int scs_sqR_Lasso_solve(double * X, double * Y, double lambda, int n, int p, double * beta);

#endif
