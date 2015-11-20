/*
 *  DESP/src/scsSRL_B.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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
#include "R_ext/BLAS.h"
#include "scsSqR_Lasso_solve.h"

SEXP scs_SRL_B(SEXP X, SEXP lambda){
	/* estimation of B */

	SEXP B, dim, status, solution, solNames;
	double * X_j, * xj, * beta;
	int i, j, k, st;
	int * stat;

	int one = 1;
	double mone = -1.0;
	int nb;
  
	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	PROTECT(solution = allocVector(VECSXP, 2));
	PROTECT(status = allocVector(INTSXP, 1));
	stat = INTEGER(status);

	/* initialize the matrix B */
	PROTECT(B = allocMatrix(REALSXP,p,p));

	X_j = (double *) Calloc((unsigned) n*(p-1), double);
	xj = (double *) Calloc((unsigned) n, double);
	beta = (double *) Calloc((unsigned) p, double);

	/* initialize status, 0 means success, 1 means failure */
	stat[0] = 0;

	/* compute an estimator of each column by the square-root Lasso */
	for(j=0; j<p; j++){
			/*for(i=0; i<n; i++){
				for(k=0; k<p; k++){
					if(k<j){
						X_j[i+n*k] = REAL(X)[i+n*k];
					}
					else if(k>j){
						X_j[i+n*(k-1)] = REAL(X)[i+n*k];
					}
				}
				xj[i] = REAL(X)[i+n*j];
			}*/
			nb = n*j;
			F77_NAME(dcopy)(&nb, REAL(X), &one, X_j, &one);
			nb = n*(p-j-1);
			F77_NAME(dcopy)(&nb, REAL(X)+(n*(j+1)), &one, X_j+(n*j), &one);
			F77_NAME(dcopy)(&n, REAL(X)+(n*j), &one, xj, &one);

			st = scs_sqR_Lasso_solve(X_j,xj,*(REAL(lambda)), n, p-1, beta);
			if (st != SCS_SOLVED){
				stat[0] = 1;
			}
			/*for(i=0; i<p; i++){
				if(i<j){
					REAL(B)[i+p*j] = - beta[i];
				}
				else if(i>j){
					REAL(B)[i+p*j] = - beta[i-1];
				}
				else REAL(B)[i+p*j] = 1;
			}*/
			F77_NAME(daxpy)(&j, &mone, beta, &one, REAL(B)+(p*j), &one);
			REAL(B)[i+p*j] = 1;
			nb = p-j-1;
			F77_NAME(daxpy)(&nb, &mone, beta, &one, REAL(B)+(p*j + j+1), &one);
	}

	SET_VECTOR_ELT(solution,0,B);
	SET_VECTOR_ELT(solution,1,status);
	
	PROTECT(solNames = allocVector(STRSXP,2));
	SET_STRING_ELT(solNames,0,mkChar("B"));
	SET_STRING_ELT(solNames,1,mkChar("status"));
	setAttrib(solution,R_NamesSymbol,solNames);

	/* free memory */
	if(xj) Free(xj);
	if(X_j) Free(X_j);
	if(beta) Free(beta);

	UNPROTECT(4);
	return solution;
}




