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

#ifndef CDSQR_LASSO
#define CDSQR_LASSO

#include <math.h>

int cd_sqR_Lasso_solve_loop(double * X, double * theta, double lambda, int n, int p, double * beta, int col, int sto);

void sample(int n, int * s);

#endif
