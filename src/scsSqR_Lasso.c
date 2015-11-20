/*
 *  DESP/src/scsSqR_Lasso.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

#include <R.h>
#include <Rinternals.h>
#include "scsSqR_Lasso_solve.h"

SEXP scs_sqR_Lasso(SEXP X, SEXP Y, SEXP lambda){

	double * bbeta;
	int * stat;
	SEXP dim, status, beta, solution, solNames;

	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	PROTECT(solution = allocVector(VECSXP, 2));
	PROTECT(beta = allocVector(REALSXP, p));
	PROTECT(status = allocVector(INTSXP, 1));
	bbeta = REAL(beta);
	stat = INTEGER(status);

	stat[0] = scs_sqR_Lasso_solve(REAL(X), REAL(Y), *(REAL(lambda)), n, p, bbeta);

	SET_VECTOR_ELT(solution,0,beta);
	SET_VECTOR_ELT(solution,1,status);
	
	PROTECT(solNames = allocVector(STRSXP,2));
	SET_STRING_ELT(solNames,0,mkChar("beta"));
	SET_STRING_ELT(solNames,1,mkChar("status"));
	setAttrib(solution,R_NamesSymbol,solNames);

	UNPROTECT(4);

	return solution;
}

